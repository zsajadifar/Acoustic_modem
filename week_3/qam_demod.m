function y = qam_demod(x , Nq)

M = 2^(Nq);
y = qamdemod(x,M, 'OutputType', 'bit','UnitAveragePower',true);