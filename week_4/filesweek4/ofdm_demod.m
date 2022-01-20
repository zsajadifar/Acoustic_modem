function y_ofdm = ofdm_demod(x_ofdm , N, cyclic_prefix , H , reminder)

x_reshape = reshape(x_ofdm ,N+cyclic_prefix ,[]);
x_reshape_without_prefix = x_reshape(cyclic_prefix+1:end,:);
y = fft(x_reshape_without_prefix ,N,1);
y1 = y ./ H; %  H is equal to  fft(h, N);
y2 = y1(2:N/2,:);
y_ofdm = y2(:);

if reminder~=0
    y_ofdm =y_ofdm(1:end-((N/2 -1)-reminder)); 
end

end
