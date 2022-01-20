function [y_ofdm_demod_on_off,H] = ofdm_demod_on_off_with_training(x_ofdm_on_off,N,cp,trainpacket,Lt,Ld,reminder_signal,p,index)

x_reshape = reshape(x_ofdm_on_off ,N+cp ,[]);
x_reshape_without_prefix = x_reshape(cp+1:end,:);
Y = fft(x_reshape_without_prefix ,N,1);

indx_signal = 1:1:size(Y,2);
indx_train=(1:1:Lt)';
% indx = 1:1:size(Y,2);
for i=1:p
    indx_train = [indx_train , indx_train(:,i) + (Lt+Ld)];
end

indx_signal(indx_train(:)) = [];
Rx_trainpacket = Y(:,indx_train(:));
Rx_signalpacket= Y(:,indx_signal);

H=zeros(N,p+1);
for i = 1 : p+1
    for j=1:N
      H(j,i) = trainpacket(j,:).'\Rx_trainpacket(j,(i-1)*Lt+1:i*Lt).';
    end
end
H(1,:) = 0;
H(N/2+1,:) = 0;


y= zeros(size(Rx_signalpacket));
for i=0:p-1
    y(:,(i*Ld+1:i*Ld+Ld))= Rx_signalpacket(:,(i*Ld+1:i*Ld+Ld))./H(:,i+1);
end
y(:,p*Ld+1:end)= Rx_signalpacket(:,p*Ld+1:end)./H(:,p+1);
y2 = y(2:N/2,:);
y2 = y2(index,:);
y_ofdm_demod_on_off = y2(:);

N_new    = length(index);
if reminder_signal~=0
    y_ofdm_demod_on_off =y_ofdm_demod_on_off(1:end-((N_new)-reminder_signal)); 
end

end