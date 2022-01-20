fs = 16000;
fsine = 1500;
t = 0: 1/fs : 2-1/fs;
sinewave = sin(2*pi*fsine*t);

plot(t(1:1000),sinewave(1:1000))

[simin,nbsecs,fs]=initparams(sinewave,fs);

sim('recplay')
out=simout.signals.values;

sound(out ,fs)
