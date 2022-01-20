close all
clear 

fs=16000;
N=1024;
noise = wgn(1,2*fs,0);
w1 = 2*700/fs;
w2 = 2*3000/fs;
b = fir1(24,[w1,w2],'stop');
figure,freqz(b)
noise_filtered = filter(b,1,noise);

[S1 , f1 , ] = spectrogram(noise,N,[],N,fs,'yaxis');
[S2 , f2 , ] = spectrogram(noise_filtered ,N,[],N,fs,'yaxis');

figure,
PSD_out_noise = mean(abs(S1).^2 , 2);
PSD_out_conv  = mean(abs(S2).^2 , 2);
subplot(2,1,1), plot(f1 ,10*log10(PSD_out_noise))
title("PSD of noise")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
ylim([-50 50])
subplot(2,1,2), plot(f2 ,20*log10(PSD_out_conv))
title("PSD  of stopbanded noise")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
ylim([-50 50])

%% IR estimation

[simin,nbsecs,fs]=initparams(noise_filtered ,fs);
sim('recplay')
out=simout.signals.values;
sound(out ,fs)

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

c = x_new ;
r = [x_new(1),zeros(imp_length-1,1)'];
x_toeplitz = toeplitz(c ,r);

h = x_toeplitz\y_new;

H= fft(h);
f = fs*(0:(length(h)/2))/length(h);
H1 = H(1:length(h)/2+1);

figure,
subplot(2,1,1), plot(h)
title("estimated IR2 with stopbanded noise")
xlabel("Sample")
ylim([-0.4 0.4])
subplot(2,1,2), plot(f,20*log10(abs(H1)))
title("frequency response of IR2 with stopbanded noise")
xlabel("Freq (Hz)")
ylabel("|H(f)| (dB)")
ylim([-40 10])

save('IRest.mat','h')
