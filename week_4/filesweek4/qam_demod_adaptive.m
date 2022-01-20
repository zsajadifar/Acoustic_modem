function y_demod = qam_demod_adaptive(x_qam , bk , reminder)

bk=bk(:);
bk_nonzero_index = find(bk~=0);
bk_no_zero = bk(bk~=0);
M = 2.^(bk_no_zero);
A = reshape(x_qam ,numel(bk),[]);

B=[];
y=[];
for j=1:numel(bk_no_zero)
    B = qamdemod(A(bk_nonzero_index(j),:),M(j), 'OutputType', 'bit','UnitAveragePower',true);
    y = [y;B];
end

y_demod = y(:);

if reminder~=0
    y_demod =y_demod(sum(bk)-reminder+1:end); 
end

end