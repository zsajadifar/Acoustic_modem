close all
clear all

%% white noise experiment 1.2
fs=16000;
N=1024;
noise = wgn(1,2*fs,0);
[simin,nbsecs,fs]=initparams(noise,fs);
sim('recplay')
out=simout.signals.values;
sound(out ,fs)

% Spectrogram 
figure,
subplot(2,1,1), spectrogram(noise,N,[],N,fs,'yaxis');
title("Spectrogram of transmitted noise")
caxis([-160 -20])
subplot(2,1,2), spectrogram(out,N,[],N,fs,'yaxis');
title("Spectrogram of recorded noise")
caxis([-160 -20])

% PSD Pwelch 
[S_noise , f_noise , ] =spectrogram(noise,N,[],N,fs,'yaxis');
[S_out   , f_out   , ] =spectrogram(out,N,[],N,fs,'yaxis');
figure,
PSD_pwelch_noise = mean(abs(S_noise).^2 , 2);
PSD_pwelch_out   = mean(abs(S_out).^2 , 2);
subplot(2,1,1), plot(f_noise ,10*log10(PSD_pwelch_noise))
title("PSD(Pwelch) of transmitted noise")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
ylim([-60 60])
subplot(2,1,2), plot(f_out ,10*log10(PSD_pwelch_out))
title("PSD(Pwelch) of recorded noise")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
ylim([-60 60])

%% IR2
x = simin(:,1);
y = out;
imp_length = 400;
fs = 16000;
delay1 = 2*fs;  % 2 Sec silence 
delay2 = fs;    % 1 Sec silence 
delay3 = 200;   % delay to make sure that the channel is casual
indx = find(y >=0.05); 
x_new  = x(delay1+1:end-delay2);
y_new = y(indx(1)-delay3+1:indx(1)-delay3+length(x_new));

% Toeplitz of x and implementing LSM
c = x_new ;
r = [x_new(1),zeros(imp_length-1,1)'];
x_toeplitz = toeplitz(c ,r);
h = x_toeplitz\y_new;

H  = fft(h);
f2 = fs*(0:(length(h)/2))/length(h);
H2 = H(1:length(h)/2+1);

figure,
subplot(2,1,1), plot(h)
title("estimated IR2")
xlabel("Sample")
ylim([-0.4 0.4])
subplot(2,1,2), plot(f2 , 20*log10(abs(H2)))
title("frequency response of IR2")
xlabel("Freq (Hz)")
ylabel("|H(f)| (dB)")
ylim([-40 10])

%% IR1
size_imp=2*fs;
imp = [1; zeros(size_imp-1,1)];
[simin,nbsecs,fs] = initparams(imp,fs);
sim('recplay')
out_imp = simout.signals.values;
sound(out_imp,fs)

fft_out_imp = fft(out_imp);
f1 = fs*(0:(length(out_imp)/2))/length(out_imp);
H1 = fft_out_imp(1:length(out_imp)/2+1);

figure,
subplot(2,1,1), plot(out_imp)
title("estimated IR1")
xlabel("Sample")
ylim([-0.4 0.4])
subplot(2,1,2), plot(f1 , 20*log10(abs(H1)))
title("frequency response of IR1")
xlabel("Freq (Hz)")
ylabel("|H(f)| (dB)")
ylim([-40 10])