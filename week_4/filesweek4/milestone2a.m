clc
close all
clear all

SNR =[5,8,10,12,15,18,20,25,30,35];
NQ  =[4 , 5 , 6 , 8] ;
string = {};

L = 1000;
x = randi([0 1], 1,L)';

for i = 1:numel(NQ)
    Nq=NQ(i);
    M = 2^(Nq);
    [x_mod,rem_qam] = qam_mod(x,Nq);
    avgPower = mean(abs(x_mod).^2);
    
    for j = 1:numel(SNR)
        
        snr         = SNR(j);
        x_mod_noise = awgn(x_mod ,snr);
        x_demod     = qam_demod(x_mod_noise , Nq ,rem_qam);
        BER(j) = ber(x , x_demod);
    end
    
    hold on
    plot(SNR , BER , 'LineWidth',2)
    string{end+1}=char(strcat('M = ',num2str(M)));
end

hold off
xlabel('SNR')
ylabel('BER')
title('BER vs SNR plot for different M-QAM')
lgd = legend(string);
lgd.NumColumns = 1;
lgd.FontWeight = 'bold';
lgd.FontSize = 12;


