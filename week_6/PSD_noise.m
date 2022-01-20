%% calculate noise PSD, ture the volume off

fs=44000;
N = 1024;
x = randn(1,2*fs);
[simin,nbsecs,fs]=initparams(x,fs);
sim('recplay')
out_noise=simout.signals.values;
% sound(out_noise ,fs)

noise = out_noise(2*fs-20+1:2*fs-20+2*fs);
[S_noise, f_noise, ]= spectrogram(noise,N,[],N,fs,'yaxis');
noise_PSD = mean(abs(S_noise).^2 , 2);
%[noise_PSD, f_noise] = pwelch(noise,N,0,N,fs);
% figure, plot(f_noise ,10*log10(noise_PSD))
figure,plot(10*log10(noise_PSD))
title("PSD noise")
xlabel("Frequency bins")
ylabel("PSD (dB)")

save('PSD_noise.mat','noise_PSD')
