function x_ofdm = ofdm_mod(x_mod , N , zero_cyclic_prefix)

P = round(length(x_mod)/((N/2)-1));
packet = zeros(N, P);

for i = 0:P-1
    Frame_first_half  = x_mod(1+((i)*(N/2 -1)):(i+1)*(N/2 -1));
    Frame_second_half = flip(conj(x_mod(1+((i)*(N/2 -1)):(i+1)*(N/2 -1))));
    packet(:,i+1) = [0;Frame_first_half;0;Frame_second_half];
end


packet_ifft           = ifft(packet,N,1);
% A                   = zeros(zero_cyclic_prefix , P);
% packet_ifft_prefix  = cat(1 , A , packet_ifft );

packet_ifft_prefix  = [packet_ifft(end-zero_cyclic_prefix+1:end,:);packet_ifft];

x_ofdm = packet_ifft_prefix(:);


