clc
clear all
close all

cp = 500 + 20; % should be: N > L > channel impulse length
N=1024;
Nq=4;
BWusage=30;
N_new = (BWusage/100)*(N/2 -1);

%% channel estimation by dummy transmission

Ld=50; % number of data frames
Lt=20; % number of training frames
num_frame=50;
x = randi([0 1], 1,Nq * (N/2 -1)); % pilot signal
x = x(:);
x_repmat = repmat(x , num_frame , 1);
[trainblock,rem_qam] = qam_mod(x,Nq);
trainblock = repmat(trainblock , num_frame , 1);
[Tx,rem_ofdm] = ofdm_mod(trainblock, N ,cp);


fs=44000;
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


%% transmite image, using estimated channel response
%% Convert BMP image to bitstream and QAM modilation
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
[qamStream ,rem_qam_signal] = qam_mod(bitStream,Nq);% symbol

%% OFDM

x = randi([0 1], 1,Nq * (N/2 -1)); % pilot signal
x = x(:);
x_repmat = repmat(x , Lt , 1);
[trainblock,rem_qam] = qam_mod(x,Nq);
trainblock = repmat(trainblock , Lt , 1);

[ofdmStream,trainpacket, rem_ofdm,p] = ofdm_mod_on_off_with_training(qamStream,N,cp,trainblock,Lt,Ld,indx);

%% channel

% Rx = simulate_channel(ofdmStream, Lt, Ld, N, cp);

fs=16000;
f=700;
% length_h=length(h);
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
[rxQamStream,H_estimate] = ofdm_demod_on_off_with_training(Rx , N, cp ,trainpacket,Lt,Ld,rem_ofdm,p,indx);
rxBitStream = qam_demod(rxQamStream , Nq , rem_qam_signal);
BER = ber(rxBitStream,bitStream);

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure,
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;

visualize_demod_on_off(H_estimate,N,fs,imageData, imageSize, bitsPerPixel,colorMap,Ld,Lt,cp,rxBitStream,Nq,indx)

%% channel estimation plot
h_estimate = ifft(H_estimate,N);
h_estimate = h_estimate(1:length_h,:);
figure,
subplot(2,1,1), plot(h_estimate)
title("estimated IR2")
xlabel("Time")
subplot(2,1,2),plot(20*log10(abs(H_estimate)))
title("frequency response estimated IR2")
xlabel("frequency bins")
ylabel("|H_estimate(f)| (dB)")





