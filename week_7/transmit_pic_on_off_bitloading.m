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
BWusage=100;
N_new = round((BWusage/100)*(N/2 -1));

%% channel estimation by dummy transmission
num_frame=50;
x = randi([0 1], 1,Nq * (N/2 -1)); % pilot signal
x = x(:);
x_repmat = repmat(x , num_frame , 1);
[trainblock,rem_qam] = qam_mod(x,Nq);
trainblock = repmat(trainblock , num_frame , 1);
[Tx,rem_ofdm] = ofdm_mod(trainblock, N ,cp);

f=700;
% length_h=length(h);
length_h=500;
pulse_sec=0.2;
t= 0:1/fs:pulse_sec-1/fs;
pulse=10*sin(2*pi*f*t);
pulse=pulse(:);
[simin,nbsecs]=initparams_withpulse(Tx,fs,pulse,length_h);
sim('recplay')
Rx_1=simout.signals.values;
% sound(Rx ,fs)
[Rx] = alignIO(Rx_1,pulse,length_h,fs);
Rx = Rx(1:length(Tx),1); 
[y_ofdm_demod,H_estimate] = ofdm_demod_estimateC(Rx , N, cp ,Tx, rem_ofdm);
y_qam_demod = qam_demod(y_ofdm_demod , Nq , rem_qam);
BER = ber(y_qam_demod,x_repmat);
figure,plot(20*log10(abs(H_estimate)))
title("frequency response estimated IR2")
xlabel("frequency bins")
ylabel("|H_estimate(f)| (dB)")

%% finding good tones, less channel atteniuation 
num_freq_bins=N_new;
H2 = H_estimate(1:length(H_estimate)/2+1);
[H2_sorted,indx1] = sort(abs(H2),'descend');
indx2 = indx1(1:num_freq_bins);
indx  = sort(indx2,'ascend');
% indx = find(H2>= 0.60*max(abs(H2)));

%% Convert BMP image to bitstream and QAM modilation
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
[qamStream ,rem_qam_signal] = qam_mod(bitStream,Nq);% symbol

%% OFDM
x = randi([0 1], 1,Nq * (N/2 -1));
x = x(:);
[trainblock,rem_qam_training] = qam_mod(x,Nq);
[ofdmStream,rem_ofdm] = ofdm_mod_on_off(qamStream,N,cp,trainblock,Lt,indx);

%% channel
[simin,nbsecs]=initparams_withpulse(ofdmStream,fs,pulse,h_length);
sim('recplay')
Rx=simout.signals.values;
[Rx] = alignIO(Rx,pulse,h_length,fs);
Rx = Rx(1:length(ofdmStream),1); 

%% ofdm demudulation
mu = 2;
[rxQamStream,H_estimate] = ofdm_demod_on_off(Rx,N,Nq,cp,trainblock,Lt,rem_ofdm,mu,indx);
rxBitStream = qam_demod(rxQamStream , Nq , rem_qam_signal);
BER = ber(rxBitStream,bitStream);

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure,
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;

visualize_demod_on_off(H_estimate,N,fs,imageData, imageSize, bitsPerPixel,colorMap,cp,rxBitStream,Nq,indx)





