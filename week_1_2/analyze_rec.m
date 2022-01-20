close all
clear all

%% 2.2
fs  = 16000;
N = 1024;
t = 0: 1/fs : 2-1/fs;
sinewave = sin(2*pi*400*t);
sig = sinewave;

%% 2.3
fs  = 16000;
N = 256;
t = 0: 1/fs : 2-1/fs;
sinewave = sin(2*pi*400*t);
sig = sinewave;

%% 2.5
fs  = 16000;
N = 1024;
t = 0: 1/fs : 2-1/fs;
sinewave = 2 + sin(2*pi*400*t);
sig = sinewave;

%% 2.6 and 2.7, effect of clipping, comment normalization in initparams function.
fs  = 16000;
N = 1024;
t = 0: 1/fs : 2-1/fs;
sinewave = 1.5*sin(2*pi*400*t);
sig = sinewave;

%% 2.8
fs  = 16000;
N = 1024;
t = 0: 1/fs : 2-1/fs;
sinewave = sin(2*pi*100*t) + sin(2*pi*200*t) +...
           sin(2*pi*500*t) + sin(2*pi*1000*t)+...
           sin(2*pi*1500*t)+ sin(2*pi*2000*t)+...
           sin(2*pi*4000*t)+ sin(2*pi*6000*t);
sig = sinewave;

%% 2.9
fs=16000;
N=1024;
sigma =1;
mu = 0;
%noise = sigma *randn(1,2*fs)+mu;
noise = wgn(1,2*fs,0);
sig = noise;
       


%% run this section for all part
[simin,nbsecs,fs]=initparams(sig,fs);
sim('recplay')
out=simout.signals.values;
sound(out ,fs)

%% Spectrogram 

figure,
subplot(2,1,1), spectrogram(sig,N,[],N,fs,'yaxis');
title("Spectrogram for input")
subplot(2,1,2), spectrogram(out,N,[],N,fs,'yaxis');
title("Spectrogram for Output")


%% PSD Pwelch 
% figure,
% [pxx,f] = pwelch(sig,N,[],N,fs);
% subplot(2,1,1) , plot(f,10*log10(pxx))
% xlabel('Frequency (Hz)')
% ylabel('PSD (dB/Hz)')
% title("PSD Pwelch")
% [pxx,f] = pwelch(out,N,[],N,fs);
% subplot(2,1,2) , plot(f,10*log10(pxx))
% xlabel('Frequency (Hz)')
% ylabel('PSD (dB/Hz)')

[Ssig , fsig , ] =spectrogram(sig,N,[],N,fs,'yaxis');
[Sout , fout , ] =spectrogram(out,N,[],N,fs,'yaxis');

figure,
Ssig_pwelch = mean(abs(Ssig).^2 , 2);
Sout_pwelch = mean(abs(Sout).^2 , 2);
subplot(2,1,1), plot(fsig ,10*log10(Ssig_pwelch))
title("Sig Pwelch PSD")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
subplot(2,1,2), plot(fout ,10*log10(Sout_pwelch))
title("Out Pwelch PSD")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")

%% PSD Bartllets
% figure,
% [pxx,f] = pwelch(sig,bartlett(N),[],N,fs);
% subplot(2,1,1) , plot(f,10*log10(pxx))
% xlabel('Frequency (Hz)')
% ylabel('PSD (dB/Hz)')
% title("PSD Bartllets")
% [pxx,f] = pwelch(out,bartlett(N),[],N,fs);
% subplot(2,1,2) , plot(f,10*log10(pxx))
% xlabel('Frequency (Hz)')
% ylabel('PSD (dB/Hz)')

[Ssig , fsig , ] =spectrogram(sig,N,0,N,fs,'yaxis');
[Sout , fout , ] =spectrogram(out,N,0,N,fs,'yaxis');

figure,
Ssig_Bartllets = mean(abs(Ssig).^2 , 2);
Sout_Bartllets = mean(abs(Sout).^2 , 2);
subplot(2,1,1), plot(fsig , 10*log10(Ssig_Bartllets))
title("Sig Bartletts PSD")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
subplot(2,1,2), plot(fout , 10*log10(Sout_Bartllets))
title("Out Bartletts PSD")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")

%% Actual PSD 

[Ssig , fsig , ] =spectrogram(sig,rectwin(length(sig)),0,length(sig),fs,'yaxis');
[Sout , fout , ] =spectrogram(out,rectwin(length(out)),0,length(out),fs,'yaxis');

figure,
Ssig_Actual = mean(abs(Ssig).^2 , 2);
Sout_Actual = mean(abs(Sout).^2 , 2);
subplot(2,1,1), plot(fsig , 10*log10(Ssig_Actual))
title("Sig Actual PSD")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
subplot(2,1,2), plot(fout , 10*log10(Sout_Actual))
title("Out Actual PSD")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")