function info = encodeTwoLevel(publicMsg, privateMsg, receiverID, level, n)
% ENCODETWOLEVEL  Two-level QR code: public message in the QR, private
% message hidden in a texture pattern that replaces the black modules.
%
%   Steps (same as the GUI in your reference):
%   1. Public QR code from publicMsg (error-correction level L/M/Q/H)
%   2. Private message -> ASCII -> GF(2^m) -> Reed-Solomon codeword (n,k)
%   3. Codeword scrambled with a secret key (random permutation)
%   4. Scrambled codeword hidden in a texture image (LSB embedding)
%   5. Black QR modules replaced by texture modules -> final QR
%   6. Texture + key stored in the folder named receiverID
%
% Requires: Communications Toolbox (gf, rsenc).  QR encoding is pure MATLAB.

if nargin < 4, level = 'M'; end
if nargin < 5, n = 14; end

m = 8;                          % bits per symbol
k = numel(privateMsg);          % message word length
if n <= k || mod(n-k, 2) ~= 0   % n and k must differ by a positive even number
    n = k + 4;
end
if n > 2^m - 1, error('n must be <= %d for m = %d', 2^m-1, m); end

%% ---------- 1. Public QR ----------
M = qr_generate(publicMsg, level);
s = 8;                                         % pixels per module
modMask = logical(kron(M, ones(s)));                   % true where black
publicQR = uint8(~modMask) * 255;

%% ---------- 2. Message encoding with Reed-Solomon ----------
asciiVal = double(privateMsg);
msgGf    = gf(asciiVal, m);                    % Galois field representation
codeGf   = rsenc(msgGf, n, k);                 % shortened RS(n,k)
codeword = double(codeGf.x);                   % 1 x n

%% ---------- 3. Scrambling ----------
key = randi(1e6);
rng(key);  perm = randperm(n);
scrambled = codeword(perm);

%% ---------- 4. Texture + hiding ----------
N = 128;
[X, Y] = meshgrid(1:N);
cross  = abs(X - N/2) < N/8 | abs(Y - N/2) < N/8;
rng(7);                                         % fixed seed -> reference pattern
ref = uint8(30 + 35*~cross + 8*randn(N));       % dark cross on dark-grey texture

bits = reshape(dec2bin(scrambled, m).', 1, []) - '0';   % 1 x (n*m)
tex  = ref;
tex(1:numel(bits)) = bitand(tex(1:numel(bits)), uint8(254)) + uint8(bits);

%% ---------- 5. Final QR: black modules -> texture ----------
[H, W] = size(modMask);
tiled   = repmat(tex, ceil(H/N), ceil(W/N));
tiled   = tiled(1:H, 1:W);
finalQR = uint8(255 * ones(H, W));
finalQR(modMask) = tiled(modMask);

%% ---------- 6. "Send" to receiver ----------
if ~exist(receiverID, 'dir'), mkdir(receiverID); end
imwrite(tex, fullfile(receiverID, 'texture.png'));
save(fullfile(receiverID, 'key.mat'), 'key', 'n', 'k', 'm');
imwrite(ref, 'reference_pattern.png');
imwrite(publicQR, 'public_qr.png');
fid = fopen('public_msg.txt','w'); fprintf(fid,'%s',publicMsg); fclose(fid);
imwrite(finalQR,  'final_qr.png');

%% ---------- Show everything ----------
fprintf('ASCII          : %s\n', num2str(asciiVal));
fprintf('Encoded (RS)   : %s\n', num2str(codeword));
fprintf('Scrambled      : %s\n', num2str(scrambled));

figure('Name', 'Encoding', 'NumberTitle', 'off');
subplot(2,2,1); imshow(ref);      title('Texture pattern');
subplot(2,2,2); imshow(tex);      title('Message-encoded texture');
subplot(2,2,3); imshow(publicQR); title('Public QR code');
subplot(2,2,4); imshow(finalQR);  title('Final QR code');

info = struct('ascii', asciiVal, 'codeword', codeword, ...
              'scrambled', scrambled, 'n', n, 'k', k);
end