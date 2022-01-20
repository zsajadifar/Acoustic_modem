% clc
clear all
% close all
warning off

fs=44000;
length_h=500;
cp = length_h + 20; % should be: N > L > channel impulse length
N=1024;
Nq=4;
gamma =10;

%% load noise PSD
load('PSD_noise.mat');

%% channel estimation by dummy transmission

Ld=50; % number of data frames
Lt=10; % number of training frames
num_frame_dummy=50;
x = randi([0 1], 1,Nq * (N/2 -1)); % pilot signal
x = x(:);
x_repmat = repmat(x , num_frame_dummy , 1);
[trainblock,rem_qam] = qam_mod(x,Nq);
trainblock = repmat(trainblock , num_frame_dummy , 1);
[Tx,rem_ofdm] = ofdm_mod(trainblock, N ,cp);

f=700;
pulse_sec=0.2;
t= 0:1/fs:pulse_sec-1/fs;
pulse=0.8*sin(2*pi*f*t);
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
figure,plot(20*log10(abs(H_estimate(1:N/2 +1))))
title("frequency response estimated IR2")
xlabel("frequency bins")
ylabel("|H_estimate(f)| (dB)")

%% compute b(k)

H2 = H_estimate(1:length(H_estimate)/2+1);
bk = floor(log2(1+ (abs(H2).^2)./(gamma* noise_PSD)));
bk(bk>14)=14;
bk(bk<1) =0;
bk = bk(2:end-1);% ignore dc and nyquist freq bin
figure,plot(bk)
title("b(k)")
xlabel("Frequency bins")
ylabel("number of bits per symbol")



%% transmite image, using estimated channel response
%% Convert BMP image to bitstream and QAM modilation
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
[qamStream ,rem_qam_signal] = qam_mod_adaptive(bitStream,bk);% symbol

%% OFDM

x = randi([0 1], 1,Nq * (N/2 -1)); % pilot signal
x = x(:);
x_repmat = repmat(x , Lt , 1);
[trainblock,rem_qam] = qam_mod(x,Nq);
trainblock = repmat(trainblock , Lt , 1);

[ofdmStream,trainpacket, rem_ofdm,p] = ofdm_mod_with_training(qamStream,N,cp,trainblock,Lt,Ld);

%% channel

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
[rxQamStream,H_estimate] = ofdm_demod_with_training(Rx , N, cp ,trainpacket,Lt,Ld,rem_ofdm,p);
rxBitStream = qam_demod_adaptive(rxQamStream , bk , rem_qam_signal);
BER = ber(rxBitStream,bitStream);

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure,
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
sgtitle(['Adaptive bitloading with BER = ',num2str(BER)])

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





