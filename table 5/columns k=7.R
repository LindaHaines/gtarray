# columns for Table 5 k=7
rm(list=ls())
library(MASS)
library(gtarray)

# setting
nc=7
for(m in 24:30)
{
if(m==24)
{gc=c(0,1,2,4,7,15,19)+1}
if(m==25)
{gc=c(0,1,2,4,10,15,22)+1}
if(m==26)
{gc=c(0,1,2,5,7,17,20)+1}  
if(m==27)
{gc=c(0,1,2,4,8,13,19)+1}
if(m==28)
{gc=c(0,1,2,4,9,15,20)+1}
if(m==29)
{gc=c(0,1,2,4,7,17,22)+1}
if(m==30)
{gc=c(0,1,2,6,18,21,23)+1} 
  
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