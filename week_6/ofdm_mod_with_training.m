function [final_ofdm,train_packet,reminder_signal,p] = ofdm_mod_with_training(x_qam ,N,cyclic_prefix,trainblock,Lt,Ld)

%% ofdm of signal

x_qam    = x_qam(:);
reminder_signal = mod(length(x_qam),(N/2 -1));
if reminder_signal~=0
    x_mod_first = x_qam(1:end-reminder_signal);
    x_mod_last  = [x_qam(end-reminder_signal+1:end); zeros((N/2 -1)-reminder_signal,1)];
    half_packet = [reshape(x_mod_first,(N/2-1),[]),x_mod_last];
else 
    half_packet = reshape(x_qam,(N/2-1),[]);
end

P             = size(half_packet , 2);
signal_packet = [zeros(1,P);half_packet;zeros(1,P);flip(conj(half_packet),1)];

%% ofdm of trainblock

trainblock    = trainblock(:);
half_packet = reshape(trainblock,(N/2-1),[]);
% reminder_train = mod(length(trainblock),(N/2 -1));
% if reminder_train~=0
%     trainblock_first = trainblock(1:end-reminder_train);
%     trainblock_last  = [trainblock(end-reminder_train+1:end); zeros((N/2 -1)-reminder_train,1)];
%     half_packet = [trainblock(trainblock_first,(N/2-1),[]),trainblock_last];
% else 
%     half_packet = reshape(trainblock,(N/2-1),[]);
% end
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

packet_ifft        = ifft(final_packet,N,1);
packet_ifft_prefix = [packet_ifft(end-cyclic_prefix+1:end,:);packet_ifft];
final_ofdm     = packet_ifft_prefix(:);


end

