clc
close all
clear all

%% transmit picture BWusage=100;
%% variables
Ld=25; % number of data frames
Lt=10; % number of training frames
cp = 500 + 20; % should be: N > L > channel impulse length
N=1024;
Nq=4;
BWusage=100;
N_new = (BWusage/100)*(N/2 -1);
fs=44000;

%% channel estimation by dummy transmission
num_frame=50;
x = randi([0 1], 1,Nq * (N/2 -1)); % pilot signal
x = x(:);
x_repmat = repmat(x , num_frame , 1);
[trainblock,rem_qam] = qam_mod(x,Nq);
trainblock = repmat(trainblock , num_frame , 1);
[Tx,rem_ofdm] = ofdm_mod(trainblock, N ,cp);

f=700;
length_h=500;
pulse_sec=0.2;
t= 0:1/fs:pulse_sec-1/fs;
pulse=0.8*sin(2*pi*f*t);
pulse=pulse(:);
[simin,nbsecs]=initparams_withpulse(Tx,fs,pulse,length_h);
sim('recplay')
Rx_1=simout.signals.values;
[Rx] = alignIO(Rx_1,pulse,length_h,fs);
Rx = Rx(1:length(Tx),1); 
[y_ofdm_demod,H_estimate_dummy] = ofdm_demod_estimateC(Rx , N, cp ,Tx, rem_ofdm);
y_qam_demod = qam_demod(y_ofdm_demod , Nq , rem_qam);
BER = ber(y_qam_demod,x_repmat);
figure,plot(20*log10(abs(H_estimate_dummy)))
title("frequency response estimated IR2")
xlabel("frequency bins")
ylabel("|H_estimate(f)| (dB)")

%% finding good tones, less channel attenuation 
num_freq_bins=N_new;
H2 = H_estimate_dummy(1:length(H_estimate_dummy)/2+1);
[H2_sorted,indx1] = sort(abs(H2),'descend');
indx2 = indx1(1:num_freq_bins);
indx  = sort(indx2,'ascend');
indx = indx(1:end-1);

%% training signal
x = randi([0 1], 1,Nq * (N/2 -1)); % pilot signal
x = x(:);
x_repmat = repmat(x , Lt , 1);
[trainblock,rem_qam] = qam_mod(x,Nq);
trainblock = repmat(trainblock , Lt , 1);

%% Convert BMP image to bitstream and QAM modilation
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
[qamStream ,rem_qam_signal] = qam_mod(bitStream,Nq);% symbol

%% OFDM
[ofdmStream,trainpacket, rem_ofdm,p] = ofdm_mod_on_off_with_training(qamStream,N,cp,trainblock,Lt,Ld,indx);

%% channel
[simin,nbsecs]=initparams_withpulse(ofdmStream,fs,pulse,length_h);
sim('recplay')
Rx=simout.signals.values;
% sound(Rx ,fs)
[Rx] = alignIO(Rx,pulse,length_h,fs);
Rx = Rx(1:length(ofdmStream),1); 

%% ofdm demudulation
[rxQamStream,H_estimate] = ofdm_demod_on_off_with_training(Rx , N, cp ,trainpacket,Lt,Ld,rem_ofdm,p,indx);
rxBitStream = qam_demod(rxQamStream , Nq , rem_qam_signal);
fprintf('BER for first method (BWusage = 100) = ')
BER = ber(rxBitStream,bitStream)

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure,
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;

visualize_demod_on_off(H_estimate,N,fs,imageData, imageSize, bitsPerPixel,colorMap,Ld,Lt,cp,rxBitStream,Nq,indx)






pause(6)
%%on-off bitloading___________________________________________________

%% variables
BWusage=50;
N_new = (BWusage/100)*(N/2 -1);

%% finding good tones, less channel attenuation 
num_freq_bins=N_new;
H2 = H_estimate_dummy(1:length(H_estimate_dummy)/2+1);
[H2_sorted,indx1] = sort(abs(H2),'descend');
indx2 = indx1(1:num_freq_bins);
indx  = sort(indx2,'ascend');

%% OFDM
[ofdmStream,trainpacket, rem_ofdm,p] = ofdm_mod_on_off_with_training(qamStream,N,cp,trainblock,Lt,Ld,indx);

%% channel
[simin,nbsecs]=initparams_withpulse(ofdmStream,fs,pulse,length_h);
sim('recplay')
Rx=simout.signals.values;
% sound(Rx ,fs)
[Rx] = alignIO(Rx,pulse,length_h,fs);
Rx = Rx(1:length(ofdmStream),1); 

%% ofdm demudulation
[rxQamStream,H_estimate] = ofdm_demod_on_off_with_training(Rx , N, cp ,trainpacket,Lt,Ld,rem_ofdm,p,indx);
rxBitStream = qam_demod(rxQamStream , Nq , rem_qam_signal);
fprintf('BER for on-off bitloading method (BWusage = 50) = ')
BER = ber(rxBitStream,bitStream)

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure,
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;

visualize_demod_on_off(H_estimate,N,fs,imageData, imageSize, bitsPerPixel,colorMap,Ld,Lt,cp,rxBitStream,Nq,indx)



















