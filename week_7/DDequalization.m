close all
clear all

Nq = 4;
L = Nq*1000;
x = randi([0 1], 1,L);
x = x(:);
[Xk,rem_qam] = qam_mod(x,Nq);
Hk = 0.01+0.01i;
Nk = fft(wgn(1000,1,-20));
snr = 100;
%Yk = awgn(Xk*Hk,snr);
%Yk = Xk*Hk + Nk;
Yk = Xk*Hk;

delta = 0.9 + 0.9i;
Wk = 1/Hk' + delta;
MU = [0.5,1,1.5,2]; % between 0 and 2
alfa = 10^-4;
Er(1)=delta;
color=['r','b','g','m'];

for j=1:4
    for i = 2:1000
        mu=MU(j);
        X_tilda= Wk(i-1)'.*Yk(i);
        desired= qam_demod(X_tilda,Nq,rem_qam);
        X_hat  = qam_mod(desired,Nq);
        Wk(i)  = Wk(i-1) +  mu*Yk(i)*conj(X_hat - Wk(i-1)'*Yk(i))./(alfa + Yk(i)'*Yk(i));
        Er(i)  = conj(Wk(i)) - 1/Hk; 
    %     if abs(Er(i))< 10^-4
    %         break
    %     end
    end
    
    hold on 
    figure(1),plot(abs(Er(1:100)),color(j)),title('error')
    legend('mu=0.5','mu=1','mu=1.5','mu=2');
    hold off
    hold on
    figure(2),plot(abs(Wk(1:100)),color(j)),title('filter coeficient')
    legend('mu=0.5','mu=1','mu=1.5','mu=2');
    hold off
%     hold on
%     plot(repmat(abs(1/Hk), [size(Wk(1:100),2),1]),'-black')
%     hold off
end
