clc
clear all
close all


%% trainblock
h_length = 500;
cp = h_length + 10;
fs=16000;
Lt = 50;
Nq = 4;
N  = 1000;
L_bitstream = Nq * (N/2 -1);
x = randi([0 1], 1,L_bitstream);
x = x(:);
[trainblock,rem_qam_training] = qam_mod(x,Nq);


%% Convert BMP image to bitstream and QAM modilation
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
[qamStream ,rem_qam_signal] = qam_mod(bitStream,Nq);% symbol

%% OFDM
[ofdmStream,rem_ofdm] = ofdm_mod_Q7(qamStream,N,cp,trainblock,Lt);

%% channel

f=700;
pulse_sec=0.2;
t= 0:1/fs:pulse_sec-1/fs;
pulse=5*sin(2*pi*f*t);
pulse=pulse(:);
[simin,nbsecs]=initparams_withpulse(ofdmStream,fs,pulse,h_length);
sim('recplay')
Rx=simout.signals.values;
[Rx] = alignIO(Rx,pulse,h_length,fs);
Rx = Rx(1:length(ofdmStream),1); 

%% ofdm demudulation
mu = 0.5;
[rxQamStream,H_estimate] = ofdm_demod_Q7(Rx,N,Nq,cp,trainblock,Lt,rem_ofdm,mu);
rxBitStream = qam_demod(rxQamStream , Nq , rem_qam_signal);
BER = ber(rxBitStream,bitStream);

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure,
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;

visualize_demod(H_estimate,N,fs,imageData, imageSize, bitsPerPixel,colorMap,cp,rxBitStream,Nq)




