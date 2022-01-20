clc
clear all
close all

N = 1024;
L = 100;
Nq=4;
Lt =10;
Ld =20;
h1 = rand(1,L);
h2 = rand(1,L);
cp = 500+10;
[a,b,H_combined] = fixed_transmitter_side_beamformer(h1,h2,N);

% a=0;
% b=1;
%% Convert BMP image to bitstream and QAM modilation
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
[qamStream ,rem_qam_signal] = qam_mod(bitStream,Nq);% symbol

%% OFDM
x = randi([0 1], 1,Nq * (N/2 -1));
x = x(:);
[trainblock,rem_qam_training] = qam_mod(x,Nq);
[ofdmStream1,ofdmStream2,rem_ofdm_signal,p] = ofdm_mod_stereo(qamStream,N,cp,trainblock,Lt,Ld,a,b);

%% channel
Rx1 = fftfilt(h1,ofdmStream1);
Rx2 = fftfilt(h2,ofdmStream2);

Rx = Rx1 + Rx2;

snr = 50;
Rx = awgn(Rx,snr);

%% ofdm demudulation
[rxQamStream,H_estimate] = ofdm_demod_stereo(Rx,N,cp,trainblock,Lt,Ld,rem_ofdm_signal,p);
rxBitStream = qam_demod(rxQamStream , Nq , rem_qam_signal);
BER = ber(rxBitStream,bitStream);

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure,
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;







