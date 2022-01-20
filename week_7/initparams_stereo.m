function [simin,nbsecs]=initparams_stereo(toplay_left,toplay_right,fs, pulse,length_h)

toplay_left = toplay_left (:); % force column structure
toplay_right = toplay_right (:);

toplay_left = toplay_left ./ max(abs(toplay_left));
toplay_right=toplay_right ./ max(abs(toplay_right));

begin_silence = zeros(2*fs,1);
end_silence   = zeros(1*fs,1);

second_zeros = zeros(length_h,1);

simin_left= [begin_silence;pulse;second_zeros;toplay_left;end_silence];
simin_right= [begin_silence;pulse;second_zeros;toplay_right;end_silence];


% simin =simin ./ max(abs(simin));

simin(:,1) = simin_left;
simin(:,2) = simin_right;

nbsecs = length(simin)/fs;