clc
close all
clear all

% Exercise session 4: DMT-OFDM transmission scheme

%% parameters
L = 100; % arbitrary channel impulse response length
snr = 47;
cyclic_prefix=550;
N = 1024;
Nq = 4;
M = 2^(Nq);

%% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

%% QAM modulation
[qamStream ,rem_qam] = qam_mod(bitStream,Nq);% symbol
% scatterplot(qamStream)
% title(sprintf('M-QAM, M=%d',M))

%% OFDM modulation
[ofdmStream,rem_ofdm] = ofdm_mod(qamStream , N , cyclic_prefix);

%% Channel
% rxOfdmStream = ofdmStream;
% rxOfdmStream = awgn(rxOfdmStream ,snr);
% h = 2*rand(L , 1) - 1;
% H = fft(h, N);
% H=H(:);
% rx_channel = fftfilt(h,rxOfdmStream);
% rxOfdmStream = rx_channel ;

% Acoustic Channel
rxOfdmStream = ofdmStream;
rxOfdmStream = awgn(rxOfdmStream ,snr);
h_IR2 = load('IRest.mat');
h_IR2 = h_IR2.h;
H = fft(h_IR2, N);
H=H(:);
rx_channel = fftfilt(h_IR2,rxOfdmStream);
rxOfdmStream = rx_channel ;


%% OFDM demodulation
rxQamStream= ofdm_demod(rxOfdmStream, N ,cyclic_prefix , H , rem_ofdm);

%% QAM demodulation
rxBitStream = qam_demod(rxQamStream , Nq , rem_qam);

%% Compute BER
BER = ber(bitStream ,rxBitStream);

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
