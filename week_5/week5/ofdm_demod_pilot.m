function [y_ofdm_on_off,H3] = ofdm_demod_pilot(x_ofdm_on_off , N, cp , Tx , indx,reminder)

x_reshape = reshape(x_ofdm_on_off ,N+cp ,[]);
x_reshape_without_prefix = x_reshape(cp+1:end,:);
Y = fft(x_reshape_without_prefix ,N,1);

Tx_reshape = reshape(Tx ,N+cp,[]);
Tx_reshape_without_prefix = Tx_reshape(cp+1:end,:);
X = fft(Tx_reshape_without_prefix,N,1);

H=[];
for i = 1:N
    H(i) = X(i,:).'\Y(i,:).';
end
H(1) = 0;
H(N/2+1) = 0;


%% H interpolation 
H_new = H(indx+1); %%256 point oversampled to 512
H2 = [0,H_new ,0, flip(conj(H_new))];
h = ifft(H2);
H3 = fft(h,N);
h2 = ifft(H3);

figure,plot(abs(H2))
figure,plot(h)
figure,plot(abs(H3))
figure,plot(h2)


y1 = Y ./ H3.'; %  H is equal to  fft(h, N);
y2 = y1(2:N/2,:);
y2 = y2(indx,:);
y_ofdm_on_off = y2(:);

N_new = length(indx);
if reminder~=0
    y_ofdm_on_off =y_ofdm_on_off(1:end-(N_new-reminder)); 
end

end