#R
#LPBNI

#input: the known lncRNA-protein interaction adjacent matrix pro_lnc_matrix
#the row and the column of  pro_lnc_matrix represent protein and lncRNA, respectively.

#output:deci,deci(i,j) means the score of protein i for lncRNA j calculated by LPBIN

B=t(pro_lnc)
ll=nrow(B)
pp=ncol(B)
W=matrix(0,pp,pp)
deci=matrix(0,pp,ll)

for (i in 1:pp)
{
gmq=matrix(rep(B[ ,i],pp),ll,pp)*B
q=apply(gmq,2,function(gmq,a) gmq/a,a=matrix(rowSums(B),ll,1))
W[i,1:pp]=(1/colSums(B))*colSums(q)        
}
for (j in 1:ll)
{
f0=B[j, ]
f1=W%*%f0
deci[1:pp,j]=f1
}
