clear all
close all

%% 1.1
fs=16000;
size_imp=1*fs;
imp = [1; zeros(size_imp-1,1)];
[simin,nbsecs,fs] = initparams(imp,fs);
sim('recplay')
out_imp = simout.signals.values;
sound(out_imp,fs)

fft_out_imp= fft(out_imp);

figure,
subplot(2,1,1), plot(out_imp)
title("estimated impulse response")
xlabel("Sample")
ylabel("")
subplot(2,1,2), plot(10*log10(abs(fft_out_imp)))
title("frequency response")
xlabel("Freq (Hz)")
ylabel("|H(f)| (dB)")

%% 1.5
fs=16000;
N=1024;
noise_length = 1; %Second
sigma =1;
mu = 0;
%noise = sigma *randn(1,2*fs)+mu;
noise = wgn(1,noise_length*fs,0);
[simin,nbsecs,fs]=initparams(noise,fs);
sim('recplay')
out_noise=simout.signals.values;
sound(out_noise ,fs)

out_conv = fftfilt(noise',out_imp);

figure,
subplot(2,1,1), spectrogram(out_noise,N,[],N,fs,'yaxis');
title("Spectrogram of transmitted noise")
caxis([-150 -30])
subplot(2,1,2), spectrogram(out_conv,N,[],N,fs,'yaxis');
title("Spectrogram of convolved noise")
caxis([-150 -30])

[S1 , f1 , ] =spectrogram(out_noise,N,[],N,fs,'yaxis');
[S2 , f2 , ] =spectrogram(out_conv ,N,[],N,fs,'yaxis');

figure,
PSD_out_noise = mean(abs(S1).^2 , 2);
PSD_out_conv  = mean(abs(S2).^2 , 2);
subplot(2,1,1), plot(f1 ,10*log10(PSD_out_noise))
title("PSD of transmitted noise")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
ylim([-50 50])
subplot(2,1,2), plot(f2 ,10*log10(PSD_out_conv))
title("PSD  of convolved noise")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
ylim([-50 50])