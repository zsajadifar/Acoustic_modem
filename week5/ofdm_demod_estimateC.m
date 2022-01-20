function [y_ofdm,H]= ofdm_demod_estimateC(x_ofdm , N, cp ,Tx , reminder)

x_reshape = reshape(x_ofdm ,N+cp ,[]);
x_reshape_without_prefix = x_reshape(cp+1:end,:);
Y  = fft(x_reshape_without_prefix ,N,1);

Tx_reshape = reshape(Tx ,N+cp,[]);
Tx_reshape_without_prefix = Tx_reshape(cp+1:end,:);
X = fft(Tx_reshape_without_prefix,N,1);

for i = 1 : N
    H(i) = X(i,:).'\Y(i,:).';
end
H(1) = 0;
H(N/2+1) = 0;
H = H(:);
y1 = Y ./ H; %H is equal to  fft(h, N);
y2 = y1(2:N/2,:);
y_ofdm = y2(:);

if reminder~=0
    y_ofdm =y_ofdm(1:end-((N/2 -1)-reminder)); 
end

end
