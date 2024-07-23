% Encode and Decode Message with Hamming Code
% Set the values of the codeword length and message length.
n = 7; % Codeword length
k = 4; % Message length
% Create a random binary message with the length equal to the message length.
data = randi([0 1], k, 1);
disp('Message:');
disp(data');
% Encode the message.
encData = encode(data, n, k, 'hamming/binary');
disp('Codeword:');
disp(encData');
% Corrupt the encoded message by introducing an error at a random location.
errLoc = randerr(1, n);
encData = mod(encData + errLoc',2);
% Decode the corrupted sequence. Observe that the decoder has correctly recovered the message.
decData = decode(encData, n, k, 'hamming/binary');
numerr = biterr(data, decData);
% Display the message, codeword, and decoded message.
disp('errCodeword:');
disp(encData');
disp('Decoded Message:');
disp(decData');

