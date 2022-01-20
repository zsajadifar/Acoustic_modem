function [ofdmStream1,ofdmStream2,reminder,p] = ofdm_mod_stereo(qamStream,N,cp,trainblock,Lt,Ld,a,b)

%% ofdm of signal

qamStream    = qamStream(:);
reminder = mod(length(qamStream),(N/2 -1));
if reminder~=0
    x_mod_first = qamStream(1:end-reminder);
    x_mod_last  = [qamStream(end-reminder+1:end); zeros((N/2 -1)-reminder,1)];
    signal_half_packet = [reshape(x_mod_first,(N/2-1),[]),x_mod_last];
else 
    signal_half_packet = reshape(qamStream,(N/2-1),[]);
end

P             = size(signal_half_packet , 2);
signal_packet = [zeros(1,P);signal_half_packet;zeros(1,P);flip(conj(signal_half_packet),1)];

%% ofdm of trainblock 

trainblock = repmat(trainblock(:) , Lt , 1);
train_half_packet = reshape(trainblock,[],Lt);

P            = size(train_half_packet , 2);
train_packet = [zeros(1,P);train_half_packet;zeros(1,P);flip(conj(train_half_packet),1)];


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

a = a(:);
b = b(:);
A = final_packet.*a;
B = final_packet.*b;
A(1,:)=0;
B(1,:)=0;
A(N/2+1,:)=0;
B(N/2+1,:)=0;

%% IFFT ,  Add cyclic prefix and parallel to serial

A_ifft        = ifft(A,N,1);
A_ifft_prefix = [A_ifft(end-cp+1:end,:);A_ifft];
ofdmStream1   = A_ifft_prefix(:);

B_ifft        = ifft(B,N,1);
B_ifft_prefix = [B_ifft(end-cp+1:end,:);B_ifft];
ofdmStream2     = B_ifft_prefix(:);


end

