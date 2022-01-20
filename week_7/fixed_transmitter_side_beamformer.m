function [a,b,H_combined] = fixed_transmitter_side_beamformer(h1,h2,N)

% H1 = fft(h1,N);
% H2 = fft(h2,N);
H1= h1;
H2= h2;
a = conj(H1)./sqrt(H1.*conj(H1)+H2.*conj(H2)); 
b = conj(H2)./sqrt(H1.*conj(H1)+H2.*conj(H2)); 
H_combined = sqrt(H1.*conj(H1)+H2.*conj(H2));

figure,
plot(20*log10(abs(H1)),'r');
hold on
plot(20*log10(abs(H2)),'b');
hold on
plot(20*log10(H_combined),'g');
legend('H left','H right','H combined')

% figure,plot(20*log10(abs(H1)));title('H left')
% figure,plot(20*log10(abs(H2)));title('H right')
% figure,plot(20*log10(H_combined));title('H combined')
end