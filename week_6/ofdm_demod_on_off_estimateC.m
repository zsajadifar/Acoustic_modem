function [y_ofdm_on_off,H] = ofdm_demod_on_off_estimateC(x_ofdm_on_off , N, cyclic_prefix , Tx , indx,reminder)

x_reshape = reshape(x_ofdm_on_off ,N+cyclic_prefix ,[]);
x_reshape_without_prefix = x_reshape(cyclic_prefix+1:end,:);
Y = fft(x_reshape_without_prefix ,N,1);

Tx_reshape = reshape(Tx ,N+cyclic_prefix,[]);
Tx_reshape_without_prefix = Tx_reshape(cyclic_prefix+1:end,:);
X = fft(Tx_reshape_without_prefix,N,1);

H=[];
for i = indx+1
    H(i) = X(i,:).'\Y(i,:).';
end
H(1) = 0;
H(N/2+1) = 0;

y1 = y ./ H; %  H is equal to  fft(h, N);
y2 = y1(2:N/2,:);
y2 = y2(indx,:);
y_ofdm_on_off = y2(:);

N_new = length(indx);
if reminder~=0
    y_ofdm_on_off =y_ofdm_on_off(1:end-(N_new-reminder)); 
end

end