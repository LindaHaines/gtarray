# columns for Table 5 k=3
rm(list=ls())
library(MASS)
library(gtarray)

# setting
nc=3
for(m in 10:20)
{
if(m==10 | m==11)
  {gc=c(0,1,3)+1}
if(m >=12 & m <= 18)
  {gc=c(0,1,4)+1}
if(m==19)
  {gc=c(0,1,8)+1}
if(m==20)
{gc=c(0,1,5)+1}

# setting
print(gc)
ntest=m*(m-nc)

# basic matrices
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