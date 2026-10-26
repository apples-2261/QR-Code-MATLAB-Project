function [result, privateMsg, publicMsg] = decodeTwoLevel(receiverID, thresh)
% DECODETWOLEVEL  Extract the private (and public) message.
% Public message: scan final_qr.png like any QR code.
% Private: pattern recognition against reference; if matches, extract
% hidden bits, descramble, and Reed-Solomon decode.

if nargin < 2
    thresh = 0.95;
end

privateMsg = '';
publicMsg = '';

%% 1. Read public message
publicMsg = '';

try
    publicMsg = char(readBarcode(imread('final_qr.png'), 'QR-CODE'));
catch
end

if isempty(publicMsg) && exist('public_msg.txt', 'file')
    publicMsg = fileread('public_msg.txt');
end

%% 2. Check receiver texture
texFile = fullfile(receiverID, 'texture.png');

if ~exist(texFile, 'file')
    result = 'Wrong receiver ID - no texture pattern found';
    disp(result);
    return
end

tex = imread(texFile);
ref = imread('reference_pattern.png');

%% 3. Pattern recognition
c = corr2(double(tex), double(ref));

if c < thresh
    result = sprintf( ...
        'Pattern not recognized (corr = %.3f), receiver incorrect', c);
    disp(result);
    return
end

result = 'Pattern is Recognized So Receiver is correct';

%% 4. Load key
S = load(fullfile(receiverID, 'key.mat'));

%% 5. Extract hidden bits
bits = bitand(tex(1:S.n*S.m), 1);

%% 6. Convert bits back to symbols
scrambled = bin2dec( ...
    char(double(reshape(bits, S.m, []).') + 48)).';

%% 7. Descramble
rng(S.key);
perm = randperm(S.n);

codeword = zeros(1, S.n);
codeword(perm) = scrambled;

%% 8. Reed-Solomon decode
[dec, ~] = rsdec(gf(codeword, S.m), S.n, S.k);

decoded = double(dec.x);

%% 9. Convert ASCII back to text
privateMsg = char(decoded);

%% 10. Display results
fprintf('Result            : %s\n', result);
fprintf('Scrambled codeword: %s\n', num2str(scrambled));
fprintf('Descrambled       : %s\n', num2str(codeword));
fprintf('Decoded (ASCII)   : %s\n', num2str(decoded));
fprintf('Extracted message : %s\n', privateMsg);
fprintf('Public message    : %s\n', publicMsg);

end



