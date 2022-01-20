function [x_ofdm,reminder] = ofdm_mod(x_mod , N , cyclic_prefix)

% size correction
x_mod    = x_mod(:);
reminder = mod(length(x_mod),(N/2 -1));
if reminder~=0
    x_mod_first = x_mod(1:end-reminder);
    x_mod_last  = [x_mod(end-reminder+1:end); zeros((N/2 -1)-reminder,1)];
    half_packet = [reshape(x_mod_first,(N/2-1),[]),x_mod_last];
else 
    half_packet = reshape(x_mod,(N/2-1),[]);
end

P = size(half_packet , 2);
packet             = [zeros(1,P);half_packet;zeros(1,P);flip(conj(half_packet),1)];
packet_ifft        = ifft(packet,N,1);
packet_ifft_prefix = [packet_ifft(end-cyclic_prefix+1:end,:);packet_ifft];
x_ofdm             = packet_ifft_prefix(:);

end

% P = round(length(x_mod)/((N/2)-1));
% packet = zeros(N, P);
% for i = 0:P-1
%     Frame_first_half  = x_mod(1+((i)*(N/2 -1)):(i+1)*(N/2 -1));
%     Frame_second_half = flip(conj(x_mod(1+((i)*(N/2 -1)):(i+1)*(N/2 -1))));
%     packet(:,i+1) = [0;Frame_first_half;0;Frame_second_half];
% end
% packet_ifft           = ifft(packet,N,1);
% % A                   = zeros(cyclic_prefix , P);
% % packet_ifft_prefix  = cat(1 , A , packet_ifft ); 
% packet_ifft_prefix    = [packet_ifft(end-cyclic_prefix+1:end,:);packet_ifft];
% 
% x_ofdm = packet_ifft_prefix(:);
% end
