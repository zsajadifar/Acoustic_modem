clc
clear all
close all


fs=44000;
length_h=500;
cp = length_h + 20; % should be: N > L > channel impulse length
N=1024;
Nq=4;
indx = 2:2:N/2 -1;
N_new=length(indx);
L = Nq * (N_new);
x = randi([0 1], 1,L);
x = x(:);
[trainblock,~] = qam_mod(x,Nq);

%% transmite image, using estimated channel response
%% Convert BMP image to bitstream and QAM modilation
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
[qamStream ,rem_qam_signal] = qam_mod(bitStream,Nq);% symbol

%% OFDM

[ofdmStream,rem_ofdm] = ofdm_mod_pilot(qamStream,N,cp,trainblock,indx);

%% channel

f=700;
pulse_sec=0.2;
t= 0:1/fs:pulse_sec-1/fs;
pulse=5*sin(2*pi*f*t);
pulse=pulse(:);
[simin,nbsecs]=initparams_withpulse(ofdmStream,fs,pulse,length_h);
sim('recplay')
Rx=simout.signals.values;
% sound(Rx ,fs)
[Rx] = alignIO(Rx,pulse,length_h,fs);
Rx = Rx(1:length(ofdmStream),1); 

%% ofdm demudulation
[rxQamStream,H_estimate] = ofdm_demod_pilot(Rx , N, cp ,trainblock,indx,rem_ofdm);
rxBitStream = qam_demod(rxQamStream , Nq , rem_qam_signal);
BER = ber(rxBitStream,bitStream);

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure,
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;

% visualize_demod_on_off(H_estimate,N,fs,imageData, imageSize, bitsPerPixel,colorMap,Ld,Lt,cp,rxBitStream,Nq,indx)
