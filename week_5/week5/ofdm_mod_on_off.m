function [x_ofdm_on_off,reminder] = ofdm_mod_on_off(x_mod , N , cp , indx)

% size correction
x_mod    = x_mod(:);
N_new    = length(indx);
reminder = mod(length(x_mod),N_new);
if reminder~=0
    x_mod_first = x_mod(1:end-reminder);
    x_mod_last  = [x_mod(end-reminder+1:end); zeros(N_new-reminder,1)];
    half_packet = [reshape(x_mod_first,N_new,[]),x_mod_last];
else 
    half_packet = reshape(x_mod,N_new,[]);
end

P = size(half_packet , 2);
half_packet_new = zeros((N/2 -1), P);
half_packet_new(indx,:) = half_packet(:,:);
% for i = 1:P
%     half_packet_new(indx,i) = half_packet(:,i) ;
% end

packet      = [zeros(1,P);half_packet_new;zeros(1,P);flip(conj(half_packet_new),1)];
packet_ifft = ifft(packet,N,1);
packet_ifft_prefix = [packet_ifft(end-cp+1:end,:);packet_ifft];
x_ofdm_on_off = packet_ifft_prefix(:);

end


