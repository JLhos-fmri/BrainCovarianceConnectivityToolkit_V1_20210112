function [R P labs] = scatcorr(A_sig,B_sig,COV)
if nargin<3
    covind = 0;
else
    covind = 1;
end
for i = 1:size(A_sig,2)
    allzero(i) = any(A_sig(:,i));
end
nonzeroind = find(allzero);
A_sig2 = A_sig(:,nonzeroind);
labst = sparse(zeros(size(A_sig2)));
labs = sparse(zeros(size(A_sig)));
Rt = zeros(size(A_sig2,2),1);
Pt = zeros(size(A_sig2,2),1);
R = zeros(size(A_sig,2),1);
P = zeros(size(A_sig,2),1);
if covind==0
    parfor i = 1:size(A_sig2,2)
        [Pi p ol] = Shepherd(A_sig2(:,i),B_sig,1000);
        labst(:,i) = sparse(ol);
        [r p] = corr(A_sig2(~ol,i),B_sig(~ol));
        Rt(i) = r;
        Pt(i) = p;
    end
    R(nonzeroind) = Rt;
    P(nonzeroind) = Pt;
    labs(:,nonzeroind) = labst;
else
    
end

end