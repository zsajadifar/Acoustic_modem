function [x_qam,reminder] = qam_mod_adaptive(x,bk)

x=x(:);
bk=bk(:);
reminder = mod(length(x),sum(bk));
if reminder~=0
    x = [zeros(sum(bk)-reminder,1);x];
end

A = reshape(x,sum(bk),[]);
bk_nonzero_index = find(bk~=0);
bk_nonzero = bk(bk~=0);
M = 2.^(bk_nonzero);
index_2 = cumsum(bk_nonzero);
index_1 = [1;index_2(1:end-1)+1];

C= zeros(length(bk),size(A,2));
for i=1:length(bk_nonzero)
    B = A(index_1(i):index_2(i),:); 
    C(bk_nonzero_index(i),:) = qammod(B(:), M(i), 'InputType', 'bit','UnitAveragePower',true);
end

x_qam = C(:);

