function [y_ofdm,H_estimate]= ofdm_demod_Q7(Rx,N,Nq,cp,trainblock,Lt,reminder,mu)

warning off 

y = reshape(Rx ,N+cp ,[]);
y2 = y(cp+1:end,:);
Yk = fft(y2(:,Lt+1:end),N,1);

%% initial estimate of channel
train_packet = y2(:,1:Lt);
Rx_train_packet  = fft(train_packet,N,1);
trainblock = repmat(trainblock , Lt , 1);
Tx_train   = reshape(trainblock,[],Lt);
Tx_train_packet= [zeros(1,Lt);Tx_train;zeros(1,Lt);flip(conj(Tx_train),1)];
H=zeros(N,1);
for j=1:N
  H(j) = Tx_train_packet(j,:).'\Rx_train_packet(j,:).';
end
H(1,:) = 0;
H(N/2+1,:) = 0;

%% NLMS 
Wk = zeros(size(Yk));
Wk(:,1)= 1./H';
Wk(1,1) = 0;
Wk(N/2+1,1) = 0;
% mu = 2; % between 0 and 2
alfa = 10^-8;

for i = 2:size(Yk,2)
    X_tilda= (Wk(:,i-1)').'.*Yk(:,i);
    desired= qamdemod(X_tilda,2^Nq, 'OutputType', 'bit','UnitAveragePower',true);
    X_hat  = qammod(desired,2^Nq, 'InputType', 'bit','UnitAveragePower',true);
    Wk(:,i)  = Wk(:,i-1) +  mu*Yk(:,i).*conj(X_hat - X_tilda)./(alfa + Yk(:,i)'*Yk(:,i));
end

X_final = (Wk').'.*Yk;
H_estimate = (1./Wk').';
H_estimate(1,:) = 0;
H_estimate(N/2+1,:) = 0;

y_ofdm = X_final(2:N/2,:);
y_ofdm = y_ofdm(:);


if reminder~=0
    y_ofdm =y_ofdm(1:end-((N/2 -1)-reminder)); 
end

end
