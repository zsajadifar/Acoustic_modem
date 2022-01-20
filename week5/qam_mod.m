function [x_qam,reminder] = qam_mod(x,Nq)

x=x(:);
reminder = mod(length(x),Nq);
if reminder~=0
    x = [zeros(Nq - reminder,1);x];
end
M = 2^(Nq);
x_qam = qammod(x,M, 'InputType', 'bit','UnitAveragePower',true);

end