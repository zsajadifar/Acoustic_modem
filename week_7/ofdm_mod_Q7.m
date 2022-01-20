function [x_ofdm,reminder] = ofdm_mod_Q7(qamStream,N,cp,trainblock,Lt)

qamStream    = qamStream(:);
reminder = mod(length(qamStream),(N/2 -1));
if reminder~=0
    x_mod_first = qamStream(1:end-reminder);
    x_mod_last  = [qamStream(end-reminder+1:end); zeros((N/2 -1)-reminder,1)];
    signal_half_packet = [reshape(x_mod_first,(N/2-1),[]),x_mod_last];
else 
    signal_half_packet = reshape(qamStream,(N/2-1),[]);
end

trainblock = repmat(trainblock , Lt , 1);
train_half_packet = reshape(trainblock,[],Lt);

half_packet = [train_half_packet,signal_half_packet];


P = size(half_packet , 2);
packet             = [zeros(1,P);half_packet;zeros(1,P);flip(conj(half_packet),1)];
packet_ifft        = ifft(packet,N,1);
packet_ifft_prefix = [packet_ifft(end-cp+1:end,:);packet_ifft];
x_ofdm             = packet_ifft_prefix(:);


end
