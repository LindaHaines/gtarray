# columns for Table 5 k=9
rm(list=ls())
library(MASS)
library(gtarray)

# setting
nc=9
for(m in 30)
{if(m==30)
{gc=c(0,1,2,3,7,11,13,16,25)+1}

# basic matrices
ntest=m*(m-nc)
xmat=grxmat(m,nc)
zmat0=grzmat0(m)
# rearrange basic zmat0
zmat=grnewzmat(m,nc,gc,zmat0)

# cmat and rank
cmat=grcmat(m,xmat,zmat)
cmatmp=ginv(cmat)

#Att
nd=nrow(cmatmp)
cmattt=cmatmp[(nc+1):nd,(nc+1):nd]
imat=diag(ntest)
jmat=matrix(1,ntest,ntest)
ijmat=imat-(jmat/ntest)
att=round(2*sum(diag(cmattt%*%ijmat))/(ntest-1),4)
print(c(m,att))
}