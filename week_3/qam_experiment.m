clc
close all
clear all

N= 256;
Nq = 4;
M = 2^(Nq);
L = 1000;
x = randi([0 1], 1,L)';
num_zeros_needed = (Nq*(N/2 -1))- mod(L , Nq*(N/2 -1));

if num_zeros_needed~=0
    t = zeros(num_zeros_needed,1);
    x2 = [t ; x];
end

x_mod = qam_mod(x2,Nq);
avgPower = mean(abs(x_mod).^2);
scatterplot(x_mod)
title(sprintf('M-QAM, M=%d',M))

snr = 20;
x_mod_noise = awgn(x_mod ,snr);
scatterplot(x_mod_noise)
title(sprintf('M-QAM, M=%d , SNR=%d ',M,snr))

x_demod = qam_demod(x_mod_noise , Nq);

x_demod_deleted_zero = x_demod(num_zeros_needed+1:end);
BER = ber(x , x_demod_deleted_zero);




