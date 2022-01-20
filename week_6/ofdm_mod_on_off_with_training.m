function [final_ofdm_on_off,train_packet,reminder_signal,p] = ofdm_mod_on_off_with_training(x_qam , N , cp ,trainblock,Lt,Ld, indx)


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
half_packet_new(indx,:) = half_packet(:,:);
% for i = 1:P
%     half_packet_new(indx,i) = half_packet(:,i) ;
% end

signal_packet = [zeros(1,P);half_packet_new;zeros(1,P);flip(conj(half_packet_new),1)];
% packet_ifft = ifft(packet,N,1);
% packet_ifft_prefix = [packet_ifft(end-cyclic_prefix+1:end,:);packet_ifft];
% x_ofdm_on_off = packet_ifft_prefix(:);


%% ofdm of trainblock 
% trainblock    = trainblock(:);
% reminder_train = mod(length(trainblock),N_new);
% if reminder_train~=0
%     trainblock_first = trainblock(1:end-reminder_train);
%     trainblock_last  = [trainblock(end-reminder_train+1:end); zeros(N_new-reminder_train,1)];
%     half_packet = [reshape(trainblock_first,N_new,[]),trainblock_last];
% else 
%     half_packet = reshape(trainblock,N_new,[]);
% end
% 
% P = size(half_packet , 2);
% half_packet_new = zeros((N/2 -1), P);
% half_packet_new(indx,:) = half_packet(:,:);
% % for i = 1:P
% %     half_packet_new(indx,i) = half_packet(:,i) ;
% % end
% 
% train_packet = [zeros(1,P);half_packet_new;zeros(1,P);flip(conj(half_packet_new),1)];

trainblock    = trainblock(:);
half_packet = reshape(trainblock,(N/2-1),[]);
P            = size(half_packet , 2);
train_packet = [zeros(1,P);half_packet;zeros(1,P);flip(conj(half_packet),1)];


%% put signal and train frames with order of Lt-Ld
num_signal_frame = size(signal_packet,2);
p = floor(num_signal_frame./Ld);
final_packet=[];

index2 = Ld:Ld:p*Ld;
index1 = [1,index2(1,1:end-1)+1];

for i=1:p
    final_packet = [final_packet,train_packet,signal_packet(: , index1(i):index2(i))];
end
final_packet = [final_packet,train_packet,signal_packet(: , index2(end)+1:end)];


%% IFFT ,  Add cyclic prefix and parallel to serial

packet_ifft          = ifft(final_packet,N,1);
packet_ifft_prefix   = [packet_ifft(end-cp+1:end,:);packet_ifft];
final_ofdm_on_off    = packet_ifft_prefix(:);


end


