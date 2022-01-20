clc
close all
clear all

cyclic_prefix=5;
snr = 20;
H = 1;
N  = 512;
Nq = 4;
M = 2^(Nq);
L = 1000;
x = randi([0 1], 1,L)';
num_zeros_needed = (Nq*(N/2 -1))- mod(L , Nq*(N/2 -1));

if num_zeros_needed~=0
    t = zeros(num_zeros_needed,1);
    x_zeropad = [t ; x];
end

x_mod = qam_mod(x_zeropad,Nq);
avgPower = mean(abs(x_mod).^2);
% scatterplot(x_mod)
% title(sprintf('M-QAM, M=%d',M))

x_ofdm_mod = ofdm_mod(x_mod , N , cyclic_prefix);


%% without noise
y_ofdm_demod = ofdm_demod(x_ofdm_mod , N ,cyclic_prefix , H);
y_ofdm_demod = qam_demod(y_ofdm_demod , Nq);

y_demod_deleted_zero = y_ofdm_demod(num_zeros_needed+1:end);
BER_without_noise = ber(x , y_demod_deleted_zero);

%% noise 

x_ofdm_mod_noise = awgn(x_ofdm_mod ,snr);
y_ofdm_demod = ofdm_demod(x_ofdm_mod_noise , N ,cyclic_prefix , H);
y_ofdm_demod = qam_demod(y_ofdm_demod , Nq);

y_demod_deleted_zero = y_ofdm_demod(num_zeros_needed+1:end);
BER_noise = ber(x , y_demod_deleted_zero);

