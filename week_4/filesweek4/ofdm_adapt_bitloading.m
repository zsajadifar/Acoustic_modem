clc
clear all
close all

%% parameters
h_IR2 = load('IRest.mat');
h_IR2 = h_IR2.h;
%snr = 47;
cyclic_prefix=length(h_IR2)+10 ;
N = 1024;
Nq = 4;
M  = 2^(Nq);
num_freq_bins = 25;
gamma = 10;

%% PSD noise 
% ture the volume off
fs = 16000;
t = 0: 1/fs : 2-1/fs;
x = randn(1,2*fs);
[simin,nbsecs,fs]=initparams(x,fs);
sim('recplay')
out_noise=simout.signals.values;
sound(out_noise ,fs)
% [S_out_noise, f_out_noise, ]= spectrogram(out_noise,N,[],N,fs,'yaxis');
% PSD_noise = mean(abs(S_out_noise).^2 , 2);
[PSD_noise, f_out_noise] = pwelch(out_noise,N,0,N,fs);
figure, 
plot(f_out_noise ,10*log10(PSD_noise))
title("PSD noise")
xlabel("Freq (Hz)")
ylabel("PSD (dB)")
ylim([-40 40])

%% channel impulse response
H = fft(h_IR2, N);
H2 = H(1:length(H)/2+1);
H=H(:);
figure,plot(h_IR2)
title("estimated IR2")
xlabel("Sample")
ylim([-0.4 0.4])
figure,plot(20*log10(abs(H2)))
title("frequency response of IR2")
xlabel("Freq (Hz)")
ylabel("|H(f)| (dB)")
ylim([-40 10])

% %find frequency bins with less channel atteniuation 
% [H2_sorted,indx1] = sort(abs(H2),'descend');
% indx2 = indx1(1:num_freq_bins);
% indx  = sort(indx2,'ascend');
indx = find(H2>= 0.80*max(abs(H2)));

%% compute b(k)
bk = floor(log2(1+ (abs(H2).^2)./(gamma* PSD_noise)));
bk(bk>14)=14;
bk(bk<1) =0;
bk = bk(2:end-1);% ignore dc and nyquist freq bin

%% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

%% Adaptive QAM modulation
[qamStream,rem_qam] = qam_mod_adaptive(bitStream,bk);

%% OFDM modulation
[ofdmStream,rem_ofdm] = ofdm_mod(qamStream , N , cyclic_prefix);
%[ofdmStream,rem_ofdm] = ofdm_mod_on_off(qamStream , N , cyclic_prefix , indx);

%% Acoustic Channel
rxOfdmStream = ofdmStream;
%rxOfdmStream = awgn(rxOfdmStream ,snr);
rx_channel = fftfilt(h_IR2,rxOfdmStream);
rxOfdmStream = rx_channel ;

%% OFDM demodulation
rxQamStream  = ofdm_demod(rxOfdmStream, N ,cyclic_prefix , H , rem_ofdm);
%rxQamStream = ofdm_demod_on_off(rxOfdmStream, N, cyclic_prefix , H , indx,rem_ofdm);

%% QAM demodulation
rxBitStream = qam_demod_adaptive(rxQamStream , bk , rem_qam );

%% Compute BER
BER = ber(bitStream ,rxBitStream);

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure,
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;




