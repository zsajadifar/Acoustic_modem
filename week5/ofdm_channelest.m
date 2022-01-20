clc
clear all
close all


h = load('IRest.mat');
h = h.h;
cp = length(h)+10; % should be: N > L > channel impulse length
N=1024;
Nq=4;
L = Nq * (N/2 -1);
x = randi([0 1], 1,L);
x = x(:);
x_repmat = repmat(x , 100 , 1);
[trainblock,rem_qam] = qam_mod(x,Nq);
pre_Tx = repmat(trainblock , 100 , 1);
[Tx,rem_ofdm] = ofdm_mod(pre_Tx, N ,cp);

Rx = fftfilt(h,Tx);
[y_ofdm_demod,H_estimate] = ofdm_demod_estimateC(Rx , N, cp ,Tx, rem_ofdm);
y_qam_demod = qam_demod(y_ofdm_demod , Nq , rem_qam);
BER = ber(y_qam_demod,x_repmat);

H = fft(h, N);
H2 = H(1:length(H)/2+1);
H=H(:);
figure,
subplot(2,1,1),plot(h)
title("IR2")
xlabel("Time")
subplot(2,1,2),plot(20*log10(abs(H)))
title("frequency response IR2")
xlabel("frequency bins")
ylabel("|H(f)| (dB)")

h_estimate = ifft(H_estimate,N);
h_estimate = h_estimate(1:length(h));
figure,
subplot(2,1,1), plot(h_estimate)
title("estimated IR2")
xlabel("Time")
subplot(2,1,2),plot(20*log10(abs(H_estimate)))
title("frequency response estimated IR2")
xlabel("frequency bins")
ylabel("|H_estimate(f)| (dB)")

%% Exercise 2.2

fs=16000;
length_h=length(h_estimate);
cp = length(h_estimate)+10; % should be: N > L > channel impulse length
N=1024;
Nq=4;
L = Nq * (N/2 -1);
x = randi([0 1], 1,L);
x = x(:);
x_repmat = repmat(x , 100 , 1);
[trainblock,rem_qam] = qam_mod(x,Nq);
pre_Tx = repmat(trainblock , 100 , 1);
[Tx,rem_ofdm] = ofdm_mod(pre_Tx, N ,cp);

f=700;
pulse_sec=0.2;
t= 0:1/fs:pulse_sec-1/fs;
pulse=10*sin(2*pi*f*t);
pulse=pulse(:);
% pulse=[ones(1,numel(t)/2),zeros(1,numel(t)/2)];
% pulse=pulse(:);
% figure,plot(t,pulse);

[simin,nbsecs]=initparams_withpulse(Tx,fs,pulse,length_h);
sim('recplay')
out=simout.signals.values;
sound(out ,fs)

[Rx] = alignIO(out,pulse,length_h,fs);
Rx = Rx(1:length(Tx),1);                %% ask
figure,plot(Rx)
title('Rx')
%Rx = [zeros(numel(Tx)-numel(Rx),1);Rx]; %% ask


[y_ofdm_demod,H_estimate] = ofdm_demod_estimateC(Rx , N, cp ,Tx, rem_ofdm);
y_qam_demod = qam_demod(y_ofdm_demod , Nq , rem_qam);
BER = ber(y_qam_demod,x_repmat);

h_estimate = ifft(H_estimate,N);
h_estimate = h_estimate(1:length_h);
figure,
subplot(2,1,1), plot(h_estimate)
title("estimated IR2")
xlabel("Time")
subplot(2,1,2),plot(20*log10(abs(H_estimate)))
title("frequency response estimated IR2")
xlabel("frequency bins")
ylabel("|H_estimate(f)| (dB)")
ylim([-20 40])




