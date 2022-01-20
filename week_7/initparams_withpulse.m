function [simin,nbsecs]=initparams_withpulse(toplay,fs, pulse,length_h)

toplay = toplay (:); % force column structure

toplay =toplay ./ max(abs(toplay));

begin_silence = zeros(2*fs,1);
end_silence   = zeros(1*fs,1);

second_zeros = zeros(length_h,1);

simin = [begin_silence;pulse;second_zeros;toplay;end_silence];

% simin =simin ./ max(abs(simin));

simin(:,2) = simin;
nbsecs = length(simin)/fs;