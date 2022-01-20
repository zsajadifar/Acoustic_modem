function [x_ofdm,reminder] = ofdm_mod_on_off(qamStream,N,cp,trainblock,Lt,indx)


%% ofdm of signal with on off loading
qamStream    = qamStream(:);
N_new    = length(indx);
reminder = mod(length(qamStream),N_new);
if reminder~=0
    x_mod_first = qamStream(1:end-reminder);
    x_mod_last  = [qamStream(end-reminder+1:end); zeros(N_new-reminder,1)];
    signal_half_packet = [reshape(x_mod_first,N_new,[]),x_mod_last];
else 
    signal_half_packet = reshape(qamStream,N_new,[]);
end

P = size(signal_half_packet , 2);
half_packet_new = zeros((N/2 -1), P);
half_packet_new(indx,:) = signal_half_packet(:,:);

trainblock = repmat(trainblock , Lt , 1);
train_half_packet = reshape(trainblock,[],Lt);
half_packet = [train_half_packet,half_packet_new];

q = size(half_packet , 2);
packet = [zeros(1,q);half_packet;zeros(1,q);flip(conj(half_packet),1)];
packet_ifft          = ifft(packet,N,1);
packet_ifft_prefix   = [packet_ifft(end-cp+1:end,:);packet_ifft];
x_ofdm    = packet_ifft_prefix(:);


end


