
C_channel = load('shannon_vs_distance.mat');
distance = 1:1:10;
figure,plot(distance , C_channel)
title("Shannon vs Distance")
xlabel("Distance (cm)")
ylabel("Channel Capacity (Kbit/Sec)")


