% wlanScrambler
% Scramble and descramble binary input sequence
% y = wlanScramble(bits,scramInit) scrambles or descrambles the binary input bits for the specified initial scramble state, using a 127-length frame-synchronous scrambler. 
% The frame-synchronous scrambler uses the generator polynomial defined in sections 17.3.5.5 and 20.3.9 of [1]. 
% The transmitter and receiver use the same scrambler to scramble bits at the transmitter and descramble bits at the receiver, respectively.
% Create the scrambler initialization and the input sequence of random bits.
scramInit = 93;
bits = randi([0,1],40,1);
% Generate the scramble sequence.
% Scrambling sequence generated using generator polynomial g(x) = x^7 + x^4 + 1.
scramInitBits = int2bit(scramInit, 7);
buffsize = min(127, size(bits, 1));
scrambleSeq = zeros(buffsize, 1, 'int8');
for i = 1:buffsize
    scrambleSeq(i) = xor(scramInitBits(1), scramInitBits(4)); % x7 xor x4
    scramInitBits(1:end-1) = scramInitBits(2:end); % Left Shift the bits
    scramInitBits(end) = scrambleSeq(i); % Update x1
end
% Replicate the scramble sequence to match the input sequence length.
scrambleSeq = repmat(scrambleSeq, ceil(size(bits, 1)/buffsize), 1);
scrambleSeq = scrambleSeq(1:size(bits, 1));
% Covert the scramble sequence to integer.
scrambleSeqInt = bit2int(scrambleSeq, 10);
% Write the scrambling sequence to a file.
fileID = fopen('scrambleSeq.txt','w');
fprintf(fileID,'%d\n', scrambleSeqInt);
fclose(fileID);
% Scramble the input sequence.
scrambleBitsTemp = xor(bits, scrambleSeq);
% Scramble the input sequence.
scrambledBits = wlanScramble(bits,scramInit);
% Test
disp(['Number of scrambled bit errors: ' num2str(sum(scrambleBitsTemp~=scrambledBits))]);
% Descramble the scrambled sequence.
descrambledBits = wlanScramble(scrambledBits, scramInit);
% Verify that the descrambled data matches the original data.
isequal(bits, descrambledBits);
% Compare the input and descrambled sequences.
disp(['Number of bit errors: ' num2str(sum(bits~=descrambledBits))]);
