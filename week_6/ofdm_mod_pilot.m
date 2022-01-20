function [x_ofdm , reminder_signal] = ofdm_mod_pilot(x_qam , N , cp ,trainblock, indx)


%% ofdm of signal with on off loading
x_qam    = x_qam(:);
N_new    = length(indx);
reminder_signal = mod(length(x_qam),N_new);
if reminder_signal~=0
    x_mod_first = x_qam(1:end-reminder_signal);
    x_mod_last  = [x_qam(end-reminder_signal+1:end); zeros(N_new-reminder_signal,1)];
    half_packet = [reshape(x_mod_first,N_new,[]),x_mod_last];
else 
    half_packet = reshape(x_qam,N_new,[]);
end

P = size(half_packet , 2);
half_packet_new = zeros((N/2 -1), P);
half_packet_new(indx,:)   = half_packet(:,:);
half_packet_new(indx-1,:) = repmat(trainblock.',P,1).';

packet = [zeros(1,P);half_packet_new;zeros(1,P);flip(conj(half_packet_new),1)];
packet_ifft = ifft(packet,N,1);
packet_ifft_prefix = [packet_ifft(end-cp+1:end,:);packet_ifft];
x_ofdm = packet_ifft_prefix(:);



end


