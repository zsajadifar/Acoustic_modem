clear all
close all

N = 1024;
fs = 16000;
K = N/2; 
B = fs/N;

%% playing no signal and record noise 
simin = 0;
nbsecs = 5;
sim('recplay')
out_noise=simout.signals.values;
sound(out_noise ,fs)

%% playing signal
pause(3)
fsine = 1500;
t = 0: 1/fs : 2-1/fs;
sig = sin(2*pi*fsine*t);
%sig = wgn(1,2*fs,0);
[simin,nbsecs,fs]=initparams(sig,fs);
sim('recplay')
out_signal_noise=simout.signals.values;
sound(out_signal_noise ,fs)

%% calculating PSD of signal 
[S_out_noise        , f_out_noise        , ]= spectrogram(out_noise,N,[],N,fs,'yaxis');
[S_out_signal_noise , f_out_signal_noise , ]= spectrogram(out_signal_noise,N,[],N,fs,'yaxis');
Ssignal = S_out_signal_noise - S_out_noise;
PSD_noise        = mean(abs(S_out_noise).^2 , 2);
PSD_signal_noise = mean(abs(S_out_signal_noise).^2 , 2);
PSD_signal       = mean(abs(Ssignal).^2 , 2);

figure, 
subplot(3,1,1), plot(f_out_noise ,10*log10(PSD_noise))
title("PSD noise")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
ylim([-40 40])
subplot(3,1,2), plot(f_out_signal_noise ,10*log10(PSD_signal_noise))
title("PSD noise+signal")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
ylim([-40 40])
subplot(3,1,3), plot(f_out_signal_noise ,10*log10(PSD_signal))
title("PSD signal")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
ylim([-40 40])

C_channel = B * sum(log2(1 + (PSD_signal./PSD_noise)));


