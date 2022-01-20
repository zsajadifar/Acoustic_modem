fs = 16000;
t = 0:1/fs :1;
imp = [1 , zeros(1,numel(t)-1)];
figure, plot((1:numel(t)),imp)
[simin,nbsecs,fs]=initparams(imp,fs);
sim('recplay')
out=simout.signals.values;
sound(out ,fs)
figure,plot((1:numel(out)),out)
fft_OUT = fft(out);
figure,plot(10*log10(abs(fft_OUT)))
title("frequency respose")
xlabel("Freq (Hz)")
ylabel("|H(f)| (dB)")


