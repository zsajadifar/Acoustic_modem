function y_ofdm = ofdm_demod(x_ofdm , N, zero_cyclic_prefix , H)

x_reshape = reshape(x_ofdm ,N+zero_cyclic_prefix ,[]);
x_reshape_without_prefix = x_reshape(zero_cyclic_prefix+1:end,:);
y = fft(x_reshape_without_prefix ,N,1);
y = y ./ H;
y1 = y(2:N/2,:);
y_ofdm = y1(:);
