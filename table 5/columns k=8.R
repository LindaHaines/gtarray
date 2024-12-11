# columns for Table 5 k=8
rm(list=ls())
library(MASS)
library(gtarray)

# setting
nc=8
for(m in 27:30)
{if(m==27)
{gc=c(0,1,2,3,6,10,16,21)+1}
if(m==28)
{gc=c(0,1,2,3,6,11,15,22)+1}
if(m==29)
{gc=c(0,1,2,4,7,12,16,22)+1}
if(m==30)
{gc=c(0,1,2,5,8,10,17,21)+1} 
  
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