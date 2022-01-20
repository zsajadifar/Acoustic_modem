function [simin,nbsecs,fs]=initparams(toplay,fs)

toplay = toplay (:); % force column structure

toplay =toplay ./ max(abs(toplay));

begin_silence = zeros(2*fs,1);
end_silence   = zeros(1*fs,1);
simin = cat(1,begin_silence,toplay,end_silence);
simin(:,2) = simin;

nbsecs = length(simin)/fs;