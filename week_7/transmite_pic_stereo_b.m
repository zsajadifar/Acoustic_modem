clc
clear all
close all

fs = 16000;
N = 1024;
L = 100;
length_h=500;
Nq=4;
Lt =10;
Ld =20;
h1 = rand(1,L);
h2 = rand(1,L);
cp = length_h+10;

%% H1 and H2 channel estimation 
x = randi([0 1], 1,Nq * (N/2 -1)); % pilot signal
x = x(:);
x_repmat = repmat(x , Lt , 1);
[trainblock,rem_qam] = qam_mod(x,Nq);
trainblock = repmat(trainblock , Lt , 1);
[Tx,rem_ofdm] = ofdm_mod(trainblock, N ,cp);

f=700;
pulse_sec=0.2;
t= 0:1/fs:pulse_sec-1/fs;
pulse=10*sin(2*pi*f*t);
pulse=pulse(:);

% pulse_sec=0.5;
% t= -pulse_sec:1/fs:pulse_sec-1/fs;
% pulse = sinc(10*t);
% pulse=pulse(:);
% figure,plot(pulse);

signal_left = [Tx ; zeros(size(Tx))];
signal_right = [zeros(size(Tx));Tx];
[simin,nbsecs]=initparams_stereo(signal_left,signal_right,fs,pulse,length_h);
sim('recplay')
Rx_1=simout.signals.values;
% sound(Rx ,fs)
[Rx] = alignIO(Rx_1,pulse,length_h,fs);
Rx = Rx(1:2*length(Tx),1); 
Rx_1 = Rx(1:length(Tx));
Rx_2 = Rx(length(Tx)+1:end);
[y_ofdm_demod_1,H1] = ofdm_demod_estimateC(Rx_1 , N, cp ,Tx, rem_ofdm);
[y_ofdm_demod_2,H2] = ofdm_demod_estimateC(Rx_2 , N, cp ,Tx, rem_ofdm);

figure,plot(20*log10(abs(H1)))
title("frequency response estimated H left")
xlabel("frequency bins")
ylabel("|H(f)| (dB)")
figure,plot(20*log10(abs(H2)))
title("frequency response estimated H rigth")
xlabel("frequency bins")
ylabel("|H(f)| (dB)")
[a,b,H_combined] = fixed_transmitter_side_beamformer(H1,H2,N);
% a=1;
% b=0;
%% Convert BMP image to bitstream and QAM modilation
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
[qamStream ,rem_qam_signal] = qam_mod(bitStream,Nq);% symbol

%% OFDM
x = randi([0 1], 1,Nq * (N/2 -1));
x = x(:);
[trainblock,rem_qam_training] = qam_mod(x,Nq);
[ofdmStream1,ofdmStream2,rem_ofdm_signal,p] = ofdm_mod_stereo(qamStream,N,cp,trainblock,Lt,Ld,a,b);

%% channel
[simin,nbsecs]=initparams_stereo(ofdmStream1,ofdmStream2,fs,pulse,length_h);
sim('recplay')
Rx_1=simout.signals.values;
% sound(Rx ,fs)
[Rx] = alignIO(Rx_1,pulse,length_h,fs);
Rx = Rx(1:length(ofdmStream1),1); 

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


