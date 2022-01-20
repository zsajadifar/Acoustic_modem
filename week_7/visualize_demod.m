function visualize_demod(H_estimate,N,fs,imageData, imageSize, bitsPerPixel,colorMap,cp,rxBitStream,Nq)


t = (N+cp)/fs;
m = (N/2 -1)*Nq;
I = [];
final_image = zeros(size(H_estimate,2),length(rxBitStream));
for k = 1:size(H_estimate,2)
    if k <size(H_estimate,2)
        I = [I;rxBitStream(1:m)];
    else
        I = [I;rxBitStream(1:end)];
    end
    final_image(k,1:length(I)) =I.';
    rxBitStream = rxBitStream(m+1:end); 
end

h_estimate = ifft(H_estimate,N);
h_estimate = h_estimate(1:500,:);
H_estimate_half = H_estimate(1:N/2+1,:);
f = fs*(0:(N/2))/N;


figure,
sec=t;
for i = 1:size(H_estimate,2) 
    imageRx = bitstreamtoimage(final_image(i,:).',imageSize, bitsPerPixel);
    subplot(2,2,1),plot(h_estimate(:,i)),title('Channel in time domain');ylim([-5 5]);drawnow;
    subplot(2,2,3),plot(f,20*log10(abs(H_estimate_half(:,i)))),title('Channel in frequency domain'),ylabel('dB'),xlabel('Frequency[Hz]');ylim([-20 30]);drawnow;
    subplot(2,2,2); colormap(colorMap); image(imageData); axis image; title('Transmitted image'); drawnow;
    subplot(2,2,4); colormap(colorMap); image(imageRx); axis image; title(['Received image after ',num2str(sec),' Sec']); drawnow;
    sec = t*(i+1);
    pause(t)
end

end
