function Q = qr_generate(str, level)
% QR_GENERATE  Pure-MATLAB QR code encoder (no toolboxes, no Java).
%   Q = qr_generate(str, level)
%     str   : text to encode (byte mode, UTF-8)
%     level : error-correction level 'L','M','Q' or 'H'  (default 'M')
%   Q : logical matrix, true = black module, includes a 4-module quiet zone.
%
%   Supports QR versions 1-10 (up to 274 bytes at level L, 122 at level H).
%   Steps: data bits -> padding -> Reed-Solomon blocks -> interleave ->
%          function patterns -> zig-zag placement -> best of 8 masks ->
%          format / version information.

if nargin < 2, level = 'M'; end
level = upper(level);
lv = find('LMQH' == level);
if isempty(lv), error('level must be L, M, Q or H'); end

% ---- tables (versions 1..10): [ecc/block, #blocks g1, data/blk g1, #blocks g2, data/blk g2]
T = [ 7 1 19 0 0; 10 1 16 0 0; 13 1 13 0 0; 17 1  9 0 0;   % v1  L M Q H
     10 1 34 0 0; 16 1 28 0 0; 22 1 22 0 0; 28 1 16 0 0;   % v2
     15 1 55 0 0; 26 1 44 0 0; 18 2 17 0 0; 22 2 13 0 0;   % v3
     20 1 80 0 0; 18 2 32 0 0; 26 2 24 0 0; 16 4  9 0 0;   % v4
     26 1 108 0 0; 24 2 43 0 0; 18 2 15 2 16; 22 2 11 2 12; % v5
     18 2 68 0 0; 16 4 27 0 0; 24 4 19 0 0; 28 4 15 0 0;   % v6
     20 2 78 0 0; 18 4 31 0 0; 18 2 14 4 15; 26 4 13 1 14; % v7
     24 2 97 0 0; 22 2 38 2 39; 22 4 18 2 19; 26 4 14 2 15; % v8
     30 2 116 0 0; 22 3 36 2 37; 20 4 16 4 17; 24 4 12 4 13; % v9
     18 2 68 2 69; 26 4 43 1 44; 24 6 19 2 20; 28 6 15 2 16]; % v10
remBitsTbl = [0 7 7 7 7 7 0 0 0 0];
alignTbl   = {[], [6 18], [6 22], [6 26], [6 30], [6 34], ...
              [6 22 38], [6 24 42], [6 26 46], [6 28 50]};
fmtLevelBits = [1 0 3 2];            % L M Q H

% ---- GF(256) tables, primitive polynomial 0x11D
expT = zeros(1, 512);  logT = zeros(1, 256);  x = 1;
for i = 0:254
    expT(i+1) = x;  logT(x+1) = i;
    x = bitshift(x, 1);
    if x >= 256, x = bitxor(x, 285); end
end
for i = 255:511, expT(i+1) = expT(i-255+1); end

% ---- choose version
data = double(unicode2native(str, 'UTF-8'));
ver = 0;
for v = 1:10
    t = T((v-1)*4 + lv, :);
    nd = t(2)*t(3) + t(4)*t(5);
    if v < 10, cc = 8; else, cc = 16; end
    if 4 + cc + 8*numel(data) <= nd*8, ver = v; break; end
end
if ver == 0, error('Text too long for versions 1-10 at level %s', level); end
t = T((ver-1)*4 + lv, :);
ecc = t(1);  nd = t(2)*t(3) + t(4)*t(5);
if ver < 10, cc = 8; else, cc = 16; end

% ---- data bit stream (byte mode)
b = [bitget(4, 4:-1:1), bitget(numel(data), cc:-1:1)];
for i = 1:numel(data), b = [b, bitget(data(i), 8:-1:1)]; end %#ok<AGROW>
b = [b, zeros(1, min(4, nd*8 - numel(b)))];          % terminator
b = [b, zeros(1, mod(-numel(b), 8))];                % byte align
cw = bin2dec(char(reshape(b, 8, []).' + 48)).';
pads = [236 17];  k = 0;
while numel(cw) < nd
    cw(end+1) = pads(mod(k,2)+1);  k = k + 1;        %#ok<AGROW>
end

% ---- Reed-Solomon per block, then interleave
lens = [repmat(t(3), 1, t(2)), repmat(t(5), 1, t(4))];
nb = numel(lens);  blocks = cell(1, nb);  ecs = cell(1, nb);  pos = 0;
g = rsGenerator(ecc, expT, logT);
for i = 1:nb
    blocks{i} = cw(pos+1 : pos+lens(i));  pos = pos + lens(i);
    ecs{i}    = rsRemainder(blocks{i}, g, expT, logT);
end
out = [];
for i = 1:max(lens)
    for j = 1:nb
        if i <= lens(j), out(end+1) = blocks{j}(i); end %#ok<AGROW>
    end
end
for i = 1:ecc
    for j = 1:nb, out(end+1) = ecs{j}(i); end %#ok<AGROW>
end
allbits = reshape(dec2bin(out, 8).' - '0', 1, []);
allbits = [allbits, zeros(1, remBitsTbl(ver))];

% ---- function patterns
sz = 17 + 4*ver;
M = zeros(sz);  F = false(sz);
corners = [0 0; 0 sz-7; sz-7 0];
for ci = 1:3
    r0 = corners(ci,1);  c0 = corners(ci,2);
    for dr = -1:7
        for dc = -1:7
            rr = r0+dr;  cc2 = c0+dc;
            if rr < 0 || cc2 < 0 || rr >= sz || cc2 >= sz, continue; end
            inside = dr>=0 && dr<=6 && dc>=0 && dc<=6;
            dk = inside && (dr==0 || dr==6 || dc==0 || dc==6 || ...
                            (dr>=2 && dr<=4 && dc>=2 && dc<=4));
            [M,F] = setm(M,F,cc2,rr,double(dk));
        end
    end
end
for i = 8:sz-9                                        % timing patterns
    [M,F] = setm(M,F,i,6,double(mod(i,2)==0));
    [M,F] = setm(M,F,6,i,double(mod(i,2)==0));
end
ap = alignTbl{ver};
for r = ap
    for c = ap
        if (r==6 && c==6) || (r==6 && c==ap(end)) || (r==ap(end) && c==6)
            continue;
        end
        for dr = -2:2
            for dc = -2:2
                [M,F] = setm(M,F,c+dc,r+dr,double(max(abs(dr),abs(dc)) ~= 1));
            end
        end
    end
end
if ver >= 7                                            % version information
    rem = ver;
    for i = 1:12, rem = bitxor(bitshift(rem,1), bitshift(rem,-11)*7973); end
    vb = bitor(bitshift(ver,12), rem);
    for i = 0:17
        bit = bitget(vb, i+1);
        a = sz - 11 + mod(i,3);  bb = floor(i/3);
        [M,F] = setm(M,F,a,bb,bit);
        [M,F] = setm(M,F,bb,a,bit);
    end
end
[M,F] = drawFormat(M,F,sz,fmtLevelBits(lv),0);         % reserve format area

% ---- zig-zag data placement
idx = 1;  up = true;  col = sz-1;
while col > 0
    if col == 6, col = col - 1; end
    if up, rows = sz-1:-1:0; else, rows = 0:sz-1; end
    for r = rows
        for c = [col, col-1]
            if ~F(r+1, c+1)
                if idx <= numel(allbits), M(r+1,c+1) = allbits(idx); end
                idx = idx + 1;
            end
        end
    end
    up = ~up;  col = col - 2;
end

% ---- choose best mask
[R, C] = ndgrid(0:sz-1, 0:sz-1);
best = inf;  bestM = M;
for mask = 0:7
    switch mask
        case 0, m = mod(R+C,2)==0;
        case 1, m = mod(R,2)==0;
        case 2, m = mod(C,3)==0;
        case 3, m = mod(R+C,3)==0;
        case 4, m = mod(floor(C/3)+floor(R/2),2)==0;
        case 5, m = mod(R.*C,2)+mod(R.*C,3)==0;
        case 6, m = mod(mod(R.*C,2)+mod(R.*C,3),2)==0;
        case 7, m = mod(mod(R+C,2)+mod(R.*C,3),2)==0;
    end
    Mt = M;
    sel = m & ~F;
    Mt(sel) = 1 - Mt(sel);
    [Mt,~] = drawFormat(Mt,F,sz,fmtLevelBits(lv),mask);   % draw into Mt!
    p = penaltyScore(Mt);
    if p < best, best = p;  bestM = Mt; end
end

% ---- add quiet zone
Q = false(sz+8);
Q(5:end-4, 5:end-4) = logical(bestM);
end

%% ------------------------------------------------------------------
function [M,F] = setm(M,F,x,y,val)   % x = column, y = row (0-based)
M(y+1, x+1) = val;  F(y+1, x+1) = true;
end

function [M,F] = drawFormat(M,F,sz,lbits,mask)
d = bitor(bitshift(lbits,3), mask);  rem = d;
for i = 1:10, rem = bitxor(bitshift(rem,1), bitshift(rem,-9)*1335); end
fb = bitxor(bitor(bitshift(d,10), rem), 21522);
gb = @(i) bitget(fb, i+1);
for i = 0:5,  [M,F] = setm(M,F,8,i,gb(i)); end
[M,F] = setm(M,F,8,7,gb(6));
[M,F] = setm(M,F,8,8,gb(7));
[M,F] = setm(M,F,7,8,gb(8));
for i = 9:14, [M,F] = setm(M,F,14-i,8,gb(i)); end
for i = 0:7,  [M,F] = setm(M,F,sz-1-i,8,gb(i)); end
for i = 8:14, [M,F] = setm(M,F,8,sz-15+i,gb(i)); end
[M,F] = setm(M,F,8,sz-8,1);                              % dark module
end

function p = gfmul(a,b,expT,logT)
if a == 0 || b == 0, p = 0; else, p = expT(logT(a+1) + logT(b+1) + 1); end
end

function g = rsGenerator(ecc, expT, logT)
g = 1;
for i = 0:ecc-1
    r = zeros(1, numel(g)+1);
    for j = 1:numel(g)
        r(j)   = bitxor(r(j),   g(j));
        r(j+1) = bitxor(r(j+1), gfmul(g(j), expT(i+1), expT, logT));
    end
    g = r;
end
end

function rem = rsRemainder(d, g, expT, logT)
ecc = numel(g) - 1;  rem = zeros(1, ecc);
for bb = d
    f = bitxor(bb, rem(1));
    rem = [rem(2:end), 0];
    for j = 1:ecc
        rem(j) = bitxor(rem(j), gfmul(g(j+1), f, expT, logT));
    end
end
end

function p = penaltyScore(M)
n = size(M,1);  p = 0;
for k = 1:n
    p = p + linePenalty(M(k,:)) + linePenalty(M(:,k).');
end
a = M(1:end-1,1:end-1);
blk = (a==M(2:end,1:end-1)) & (a==M(1:end-1,2:end)) & (a==M(2:end,2:end));
p = p + 3*sum(blk(:));
tot = n*n;  dark = sum(M(:));
p = p + 10*(ceil(abs(dark*20 - tot*10)/tot) - 1);
end

function p = linePenalty(v)
d = [1, double(diff(v) ~= 0), 1];
runs = diff(find(d));
runs = runs(runs >= 5);
p = sum(3 + runs - 5);
s = char(v + 48);
p = p + 40*(numel(strfind(s,'10111010000')) + numel(strfind(s,'00001011101')));
end

