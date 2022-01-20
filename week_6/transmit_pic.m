clc
clear all
close all
warning off

h = load('IRest.mat');
h = h.h;
cp = length(h)+20; % should be: N > L > channel impulse length
N=1024;
Nq=4;


%% Convert BMP image to bitstream and QAM modilation
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
[qamStream ,rem_qam_signal] = qam_mod(bitStream,Nq);% symbol

%% pilot signal

Ld=40; % number of data frames
Lt=50; % number of training frames
L_bitstream = Nq * (N/2 -1);
x = randi([0 1], 1,L_bitstream);
x = x(:);
[trainblock,rem_qam_training] = qam_mod(x,Nq);
trainblock = repmat(trainblock , Lt , 1);

%% OFDM
[ofdmStream,trainpacket, rem_ofdm,p] = ofdm_mod_with_training(qamStream,N,cp,trainblock,Lt,Ld);

%% channel

% Rx = simulate_channel(ofdmStream, Lt, Ld, N, cp);

fs=16000;
f=700;
length_h=length(h);
pulse_sec=0.2;
t= 0:1/fs:pulse_sec-1/fs;
pulse=10*sin(2*pi*f*t);
pulse=pulse(:);
[simin,nbsecs]=initparams_withpulse(ofdmStream,fs,pulse,length_h);
sim('recplay')
Rx=simout.signals.values;
% sound(Rx ,fs)
[Rx] = alignIO(Rx,pulse,length_h,fs);
Rx = Rx(1:length(ofdmStream),1); 

%% ofdm demudulation
[rxQamStream,H_estimate] = ofdm_demod_with_training(Rx , N, cp ,trainpacket,Lt,Ld,rem_ofdm,p);
rxBitStream = qam_demod(rxQamStream , Nq , rem_qam_signal);
BER = ber(rxBitStream,bitStream);

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure,
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;

visualize_demod(H_estimate,N,fs,imageData, imageSize, bitsPerPixel,colorMap,Ld,Lt,cp,rxBitStream,Nq)

%% channel estimation plot
h_estimate = ifft(H_estimate,N);
h_estimate = h_estimate(1:length(h),:);
figure,
subplot(2,1,1), plot(h_estimate)
title("estimated IR2")
xlabel("Time")
subplot(2,1,2),plot(20*log10(abs(H_estimate)))
title("frequency response estimated IR2")
xlabel("frequency bins")
ylabel("|H_estimate(f)| (dB)")





