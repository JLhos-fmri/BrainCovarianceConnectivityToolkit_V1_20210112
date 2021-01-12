function Zval = PtoZ(P)
disval = -65:0.001:65;
pval = normcdf(disval);
ind = find(pval<P);
Zval = disval(max(ind));
end