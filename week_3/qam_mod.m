function y = qam_mod(x,Nq)
% L = 1200;
% x = randi([0 1], 1,L)';
% Nq = 6;
M = 2^(Nq);
y = qammod(x,M, 'InputType', 'bit','UnitAveragePower',true);
% scatterplot(y)