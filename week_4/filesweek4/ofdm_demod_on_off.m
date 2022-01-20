function y_ofdm_on_off = ofdm_demod_on_off(x_ofdm_on_off , N, cyclic_prefix , H , indx,reminder)

x_reshape = reshape(x_ofdm_on_off ,N+cyclic_prefix ,[]);
x_reshape_without_prefix = x_reshape(cyclic_prefix+1:end,:);
y = fft(x_reshape_without_prefix ,N,1);
y1 = y ./ H; %  H is equal to  fft(h, N);
y2 = y1(2:N/2,:);
y2 = y2(indx,:);
y_ofdm_on_off = y2(:);

N_new = length(indx);
if reminder~=0
    y_ofdm_on_off =y_ofdm_on_off(1:end-(N_new-reminder)); 
end

end