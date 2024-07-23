% Input an image file and convert to binary stream
fileTx = 'peppers.png';
% Read image file
fData = imread(fileTx);
scale = 0.05;
origSize = size(fData);                          % Original input image size
scaledSize = max(floor(scale.*origSize(1:2)),1); % Calculate new image size
heightIx = min(round(((1:scaledSize(1))-0.5)./scale+0.5),origSize(1));
widthIx = min(round(((1:scaledSize(2))-0.5)./scale+0.5),origSize(2));
fData = fData(heightIx,widthIx,:);               % Resize image
imsize = size(fData);                            % Store new image size
txImg = fData(:);

% Plot the transmitted image
imFig.Visible = 'on';
subplot(2,1,1);
imshow(fData);
title('Transmitted Image');
subplot(2,1,2);
title('Received Image appears here...');
set(gca, 'Visible', 'off');
set(findall(gcf,'type','text'),'visible','on');

% Fragment Transmit Data
% MSDU in bytes
msduLength = 576;
numMSDUs = ceil(size(txImg, 1)/msduLength);
padZeros = msduLength - mod(size(txImg, 1), msduLength);
txData = [txImg; zeros(padZeros, 1)];
txDataBits = double(int2bit(txData, 8, false));
writematrix(txData, 'txData.txt', 'Delimiter', 'tab');

% Divide input data stream into multiple MSDUs
bitsPerOctet = 8;
data = zeros(0, 1);

for i = 0:numMSDUs-1
    % Extract image data (in Octets) for the current MSDU
    frameBody = txData(i*msduLength+1:(i+1)*msduLength, :);

    % Create MAC frame configuration object and configure sequence number
    cfgMAC = wlanMACFrameConfig(FrameType='Data',SequenceNumber=i);

    % Generate MPDU
    [psdu, lengthMPDU]= wlanMACFrame(frameBody,cfgMAC,OutputFormat='bits');
    macFrame = wlanMACFrame(frameBody, cfgMAC);
    
    % Add the MPDU to the data stream
    data = [data; psdu]; %#ok<AGROW>
end
writematrix(data, 'data.txt', 'Delimiter', 'tab');

% Generate 802.11a Baseband WLAN Signal
% Synthesize a non-HT waveform, with a non-HT format configuration object created by the wlanNonHTConfig object
cfgNonHT = wlanNonHTConfig; % Create a wlanNonHTConfig object
cfgNonHT.PSDULength = lengthMPDU; % Set the PSDU length
cfgNonHT.MCS = 6; % Set the modulation and coding scheme
cfgNonHT.ChannelBandwidth = 'CBW20'; % Set the channel bandwidth
cfgNonHT.NumTransmitAntennas = 1; % Set the number of transmit antennas
chanBW = cfgNonHT.ChannelBandwidth; % Get the channel bandwidth

% Initialize the scrambler with a random integer for each packet
scramInit = randi([1 127], numMSDUs, 1);
scramInit = [57; 65; 33]; % for test

% Set the oversampling factor to 1.5 to generate the waveform at 30 MHz for transmission
osf = 1.0;

sampleRate = wlanSampleRate(cfgNonHT); % Get the baseband sample rate

% Generate baseband NonHT packets separated by idle time
idleTime = 20e-6; % 20 microseconds
% Define default scrambler initialization value
defaultScramblerInitialization = 93;
% Windowing transition time
winTransmitTime = 1e-7; % 100 nanoseconds

useParams = struct('NumPackets', numMSDUs, 'IdleTime', idleTime, 'WindowTransitionTime', winTransmitTime, ...
    'ScramblerInitialization', scramInit, 'OversamplingFactor', osf);
overrideObjectScramInit = false;
% Validate the format configuration object is a valid type
s = validateConfig(cfgNonHT);

osf = double(useParams.OversamplingFactor); % Force from numeric to double
numUsers = 1; % Number of users
% Data packet
% Columnize and expand to a [1 Nu] cell
dataCell = repmat({int8(data(:))}, 1, numUsers);

% Number of bits in a PSDU for a single packet (convert bytes to bits)
psduLength = s.PSDULength;
numPSDUBits = psduLength*8;
% Repeat to provide initial state(s) for all users and packets
pktScramInit = scramInit(mod((0:useParams.NumPackets-1).', size(scramInit, 1))+1, :);
% Get the sampling rate of the waveform
sr = wlan.internal.cbwStr2Num(chanBW)*1e6*osf;
numTxAnt = cfgNonHT.NumTransmitAntennas;
numPktSamples = real(s.NumPPDUSamples)*osf;
% Generate the legacy preamble for NonHT packets
% Generate L-STF
%   LSTF = lstfSequence returns the sequence used for 20 MHz L-STF as 
%   per equation 20-8 in IEEE Std 802.11-2012 pg. 1695.
% Map subcarriers and replicate over subchannels
LTFS = [zeros(6,1); wlan.internal.lstfSequence(); zeros(5,1)];
% Specify FFT parameters directly for codegen
% cbw2nfft Get FFT length for the given channel bandwidth
[fftLen,numSubchannels] = wlan.internal.cbw2nfft(cfgNonHT.ChannelBandwidth);
numTones = 12*numSubchannels;
sym = repmat(LTFS, numSubchannels, 1); % Replicate for each BW
% Apply gamma rotation, replicate over antennas and apply cyclic shifts
[lstf, scalingFactor] = wlan.internal.legacyFieldMap(sym, numTones, cfgNonHT);
% OFDM modulation
modOut = wlan.internal.wlanOFDMModulate(lstf,0,osf); % 0 CP length
% Scale the preamble
lstfSym = [modOut; modOut; modOut(1:fftLen*osf/2,:)]*scalingFactor;

% Generate L-LTF 
%lltfSequence L-LTF upper and lower subcarrier sequence
%   The sequence used for 20 MHz L-LTF: equation 20-11 in IEEE Std
%   802.11-2012.
[lltfLower,lltfUpper] = wlan.internal.lltfSequence();
LLTF = [zeros(6,1); lltfLower; 0; lltfUpper; zeros(5,1)];
% Specify FFT parameters directly for codegen
[fftLen,numSubchannels] = wlan.internal.cbw2nfft(cfgNonHT.ChannelBandwidth);
numTones = 52*numSubchannels;
% Replicate L-LTF sequence for each 20MHz subchannel
sym = repmat(LLTF, numSubchannels, 1);
% Apply gamma rotation, replicate over antennas and apply cyclic shifts
[lltf, scalingFactor] = wlan.internal.legacyFieldMap(sym, numTones, cfgNonHT);
% OFDM modulate with double length GI for first symbol, no GI for second
modOut = wlan.internal.wlanOFDMModulate([lltf lltf],[fftLen/2 0],osf);
% Scale the preamble
lltfSym = modOut*scalingFactor;

% Generate L-SIG
%   Generate the non-HT SIGNAL field and non-HT OFDM transmission formats.
% Set the Rate bits
Rbits = wlan.internal.nonHTRateSignalBits(cfgNonHT.MCS);
length = cfgNonHT.PSDULength;
% Construct the SIGNAL field. Length parameter with LSB first, which is 12 bits
lengthBits = int2bit(length, 12, false);
% Even parity bit calculation
parityBit = mod(sum([Rbits; lengthBits], 1), 2);
% The SIGNAL field (IEEE Std 802.11-2016, Section 17.3.4.2)
sigBits = [Rbits; 0; lengthBits; parityBit; zeros(6, 1, 'int8')];
% Process L-SIG bits
%   Y = wlanBCCEncode(X,RATE) convolutionally encodes the binary input
%   data X using a binary convolutional code (BCC) at the specified RATE.
encodedBits = wlanBCCEncode(sigBits, '1/2');
% Interleave the encoded bits
interleavedBits = wlanBCCInterleave(encodedBits, 'Non-HT', 48);

% const = comm.internal.qam.getConstellation(Mnew, unitAveragePower);
% Use BPSK modulation for the L-SIG field
const = [-1; 1];
modData = complex(const(interleavedBits+1));
% Add pilot symbols, from IEEE Std 802.11-2016, Equation 19-14
Nsym = 1; % One symbol
z = 0;    % No offset as first symbol is with pilots
modPilot = wlan.internal.nonHTPilots(Nsym, z);
% Generate the L-SIG field
% Map subcarriers and replicate over bandwidth
ofdm = wlan.internal.vhtOFDMInfo('L-SIG', cfgNonHT.ChannelBandwidth, 1);
sym = complex(zeros(ofdm.FFTLength, 1));
sym(ofdm.ActiveFFTIndices(ofdm.DataIndices)) = repmat(modData, ofdm.NumSubchannels, 1);
sym(ofdm.ActiveFFTIndices(ofdm.PilotIndices)) = repmat(modPilot, ofdm.NumSubchannels, 1);
% Apply gamma rotation, replicate over antennas and apply cyclic shifts
[lsig, scalingFactor] = wlan.internal.legacyFieldMap(sym, ofdm.NumTones, cfgNonHT);
% OFDM modulate
cplen = ofdm.CPLength; % Cyclic prefix length
[fftLen,numSym,numTx] = size(lsig);
if osf>1
    wlan.internal.validateOFDMOSF(osf,fftLen,cplen);
    numSamples = (osf-1)*fftLen/2;        
    padding = zeros(numSamples,numSym,numTx,'like',lsig);
    lsig = [padding; lsig*osf; padding]; % Scale by OSF to account for larger FFT
    cplen = cplen*osf;
    fftLen = fftLen*osf;
end
% Shift and IFFT
postShift = ifftshift(lsig,1);
postIFFT = osf*ifft(postShift,[],1);
% Append cyclic prefix
postCP = postIFFT([end-cplen+(1:cplen),1:end],:,:,:);
symData = reshape(postCP, [(fftLen+cplen)*numSym numTx]);
% Scale the preamble
lsigSym = symData*scalingFactor;

% Assemble the preamble
preamble = [lstfSym; lltfSym; lsigSym];
psps = wlan.internal.nonhtPacketSamplesPerSymbol(cfgNonHT,s.NumDataSymbols,osf);
% Define a matrix of total simulation length
numIdleSamples = round(sr*useParams.IdleTime);
pktWithIdleLength = numPktSamples+numIdleSamples;
txWaveform = complex(zeros(useParams.NumPackets*pktWithIdleLength,numTxAnt));

% Generate the waveform for each packet
for i = 1:useParams.NumPackets
    % Extract PSDU for the current packet
    psdu = repmat({int8(1)},1,numUsers); % Cannot use cell(1, numUsers) for codegen
    for u = 1:numUsers
        psdu{u} = wlan.internal.parseInputBits(dataCell{u},numPSDUBits(u),(i-1)*numPSDUBits(u));
    end
    % Determine number of symbols and pad length
    numSym = s.NumDataSymbols;
    numPad = s.NumPadBits;
    mcsTable = wlan.internal.getRateTable(cfgNonHT);
    Nservice = 16;
    Ntail = 6;
    % Scramble padded data
    % [service; psdu; tail; pad] processing
    scramInitBits = int2bit(pktScramInit(i,:), 7);
    paddedData = [zeros(Nservice,1,'int8'); psdu{1}; zeros(Ntail,1,'int8'); zeros(numPad, 1,'int8')];
    writematrix(paddedData, 'paddedData.txt', 'Delimiter', 'tab');
    scrambData = wlanScramble(paddedData, scramInitBits);
    % Zero-out the tail bits again for encoding
    scrambData(16+size(psdu{1},1) + (1:Ntail)) = zeros(Ntail,1);
    % BCC Encoding
    encodedData = wlanBCCEncode(scrambData, mcsTable.Rate);
    % BCC Interleaving
    interleavedData = wlanBCCInterleave(encodedData, 'Non-HT', mcsTable.NCBPS);
    % Constellation mapping
    % Get symbol mapping
    symMap = comm.internal.qam.wlanSymbolMap(mcsTable.NBPSCS); 
    M = 2^mcsTable.NBPSCS;
    Mnew = double(M);   
    unitAveragePower = true;
    const = comm.internal.qam.getSquareConstellation(Mnew);
    if unitAveragePower
        defaultAveragePower = comm.internal.qam.getAveragePower(Mnew);
        const = const .* sqrt(1/defaultAveragePower);
    end
    % if input is a vector, ensure that const has same orientation as
    % input. Note that, const is a row vector.
    if iscolumn(interleavedData)
        newConst = const(:);
    else
        newConst = const;
    end
    sym = comm.internal.utilities.convertBit2Int(interleavedData, log2(Mnew));
    symbolOrderMap = coder.nullcopy(zeros(Mnew, 1, 'like', interleavedData));
    symbolOrderMap(symMap + 1) = 0:Mnew-1;
    msg = symbolOrderMap(sym + 1);
    mappedData = newConst(msg + 1);
    % Non-HT pilots, from IEEE Std 802.11-2012, Eqn 18-22
    % Reshape to form OFDM symbols
    mappedData = reshape(mappedData, mcsTable.NCBPS/mcsTable.NBPSCS, numSym);
    % Non-HT pilots, from IEEE Std 802.11-2012, Eqn 18-22
    z = 1; % Offset by 1 to account for HT-SIG pilot symbol
    pilotValues = wlan.internal.nonHTPilots(numSym,z);
    % Data packing with pilot insertion, replicate over 20 MHz subchannels
    ofdm = wlan.internal.vhtOFDMInfo('NonHT-Data', cfgNonHT.ChannelBandwidth, 1);
    packedData = complex(zeros(ofdm.FFTLength, numSym));
    packedData(ofdm.ActiveFFTIndices(ofdm.DataIndices), :) = repmat(mappedData, ofdm.NumSubchannels, 1);
    packedData(ofdm.ActiveFFTIndices(ofdm.PilotIndices), :) = repmat(pilotValues, ofdm.NumSubchannels, 1);
    % Apply gamma rotation, replicate over antennas and apply cyclic shifts
    [data,scalingFactor] = wlan.internal.legacyFieldMap(packedData, ofdm.NumTones, cfgNonHT);
    % OFDM modulate
    cplen = ofdm.CPLength; % Cyclic prefix length
    [fftLen,numSym,numTx] = size(data);
    if osf>1
        wlan.internal.validateOFDMOSF(osf,fftLen,cplen);
        numSamples = (osf-1)*fftLen/2;        
        padding = zeros(numSamples,numSym,numTx,'like',data);
        data = [padding; data*osf; padding]; % Scale by OSF to account for larger FFT
        cplen = cplen*osf;
        fftLen = fftLen*osf;
    end
    % Shift and IFFT
    postShift = ifftshift(data,1);
    postIFFT = osf*ifft(postShift,[],1); 
    % Append cyclic prefix
    postCP = postIFFT([end-cplen+(1:cplen),1:end],:,:,:);
    symData = reshape(postCP, [(fftLen+cplen)*numSym numTx]);
    % Scale the data
    dataSym = symData*scalingFactor;

    % Construct packet from preamble and data
    pktIdx = (i-1)*pktWithIdleLength+(1:numPktSamples);
    txWaveform(pktIdx,:) = [preamble; dataSym];
end

% Window waveform
wLength = 2*ceil(useParams.WindowTransitionTime*sr/2); % Window length in samples
% Ensure window transition time is within 2x smallest CP length
maxWLength = 2*min(psps.CPPerSymbol(psps.CPPerSymbol>0));
txWaveform = wlan.internal.windowWaveform(txWaveform,psps.NumSamplesPerSymbol,psps.CPPerSymbol,wLength,useParams.NumPackets,numIdleSamples);
writematrix(txWaveform, 'txWaveform.txt', 'Delimiter', 'tab');
