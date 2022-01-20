fs=16000;
noise = wgn(1,2*fs,0);

[simin,nbsecs,fs]=initparams(noise,fs);
sim('recplay')
out=simout.signals.values;
sound(out ,fs)
figure, plot(out)

x = simin(:,1);
y = out;
imp_length = 500;
fs = 16000;
delay1 = 2*fs;  % 2 Sec silence sample
delay2 = fs;    % 1 Sec silence sample
x_new  = x(delay1+1:end-delay2);

indx = find(y >=0.25*max(y));
y_new = y(indx(1)-200+1:indx(1)-200+length(x_new));


c = x_new ;
r = [x_new(1),zeros(imp_length-1,1)'];
x_toeplitz = toeplitz(c ,r);

h = x_toeplitz\y_new;

H= fft(h);

figure,
subplot(2,1,1), plot(h)
title("estimated impulse response")
xlabel("Sample")
ylabel("")
subplot(2,1,2), plot(10*log10(abs(H)))
title("frequency response")
xlabel("Freq (Hz)")
ylabel("|H(f)| (dB)")

save('IRest.mat','h')


    


