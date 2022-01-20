function [y_ofdm_on_off,H3] = ofdm_demod_pilot(Rx , N, cp , trainblock , indx,reminder)

x_reshape = reshape(Rx ,N+cp ,[]);
x_reshape_without_prefix = x_reshape(cp+1:end,:);
Y = fft(x_reshape_without_prefix ,N,1);
P = size(Y,2);

X = repmat(trainblock.',P,1).';
M = numel(indx);

H=[];
for i=1:M
    H(i) = X(M,:).'\Y(indx(i)-1,:).';
end
H(1) = 0;


%% H interpolation 
H2 = [0,H ,0, flip(conj(H))];
h = ifft(H2);
H3 = fft(h,N);
h2 = ifft(H3);

% figure,plot(abs(H2))
% figure,plot(h)
% figure,plot(abs(H3))
% figure,plot(h2)


y1 = Y ./ H3.'; %  H is equal to  fft(h, N);
y2 = y1(2:N/2,:);
y2 = y2(indx,:);
y_ofdm_on_off = y2(:);

N_new = length(indx);
if reminder~=0
    y_ofdm_on_off =y_ofdm_on_off(1:end-(N_new-reminder)); 
end

end