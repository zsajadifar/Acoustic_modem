clc
clear all
close all

fs=16000;
length_h=500;
cp =length_h +10; % should be: N > L > channel impulse length
N=1024;
Nq=4;
indx = 2:2:N/2 -1;
N_new=length(indx);
L = Nq * (N_new);
x = randi([0 1], 1,L);
x = x(:);
x_repmat = repmat(x , 100 , 1);
[trainblock,rem_qam] = qam_mod(x,Nq);
pre_Tx = repmat(trainblock , 100 , 1);
[Tx,rem_ofdm] = ofdm_mod_on_off(pre_Tx,N ,cp, indx);

f=700;
pulse_sec=0.2;
t= 0:1/fs:pulse_sec-1/fs;
pulse=0.8*sin(2*pi*f*t);
pulse=pulse(:);

[simin,nbsecs]=initparams_withpulse(Tx,fs,pulse,length_h);
sim('recplay')
out=simout.signals.values;
sound(out ,fs)

[Rx] = alignIO(out,pulse,length_h,fs);
Rx = Rx(1:length(Tx),1);                %% ask
figure,plot(Rx)
title('Rx')
%Rx = [zeros(numel(Tx)-numel(Rx),1);Rx]; %% ask



[y_ofdm_demod,H_estimate] = ofdm_demod_pilot(Rx , N, cp ,Tx, indx, rem_ofdm);
y_qam_demod = qam_demod(y_ofdm_demod , Nq , rem_qam);
BER = ber(y_qam_demod,x_repmat);

h_estimate = ifft(H_estimate,N);
h_estimate = h_estimate(1:length_h);
figure,
subplot(2,1,1), plot(h_estimate)
title("estimated IR2")
xlabel("Time")
subplot(2,1,2),plot((abs(H_estimate)))
title("frequency response estimated IR2")
xlabel("frequency bins")
ylabel("|H_estimate(f)|")
