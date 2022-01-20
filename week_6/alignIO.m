function [out_aligned] = alignIO(out,pulse,length_h,fs)

[r,lag] = xcorr(out,pulse);
[max_corr,indx] = max(abs(r));
delay = lag(indx);
out_aligned = out(delay+length(pulse)+length_h-60+1:end);
figure,plot(out)
hold on
plot(delay,out(delay),'r*')
title('Out signal with the point calculated from cross crrolation as delay')


% test=simin(:,1);
% [r,lag] = xcorr(test,pulse);
% [max_corr,indx] = max(abs(r));
% delay = lag(indx);
% out_aligned = test(delay+length(pulse)+length_h-20+1:end-fs);
% figure,plot(test)
% hold on
% plot(delay,test(delay),'r*')
% hold on 
% plot(out_aligned,'r')
