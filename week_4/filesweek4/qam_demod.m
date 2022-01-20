function y = qam_demod(x_qam , Nq , reminder)

M = 2^(Nq);
y = qamdemod(x_qam,M, 'OutputType', 'bit','UnitAveragePower',true);

if reminder~=0
    y =y(Nq-reminder+1:end); 
end

end