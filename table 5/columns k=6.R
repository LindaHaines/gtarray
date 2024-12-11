# columns for Table 5 k=6
rm(list=ls())
library(MASS)
library(gtarray)

# setting
nc=6
for(m in 20:30)
{if(m==20 | m==21)
{gc=c(0,1,2,4,10,15)+1}
if(m==22)
{gc=c(0,1,2,6,9,11)+1}
if(m==23)
{gc=c(0,1,2,5,11,16)+1}
if(m==24)
{gc=c(0,1,2,6,9,14)+1}
if(m==25)
{gc=c(0,1,2,6,9,15)+1}
if(m==26)
{gc=c(0,1,2,5,12,18)+1}
if(m==27)
{gc=c(0,1,2,5,13,22)+1}  
if(m==28)
{gc=c(0,1,4,15,20,22)+1}
if(m==29)
{gc=c(0,1,2,6,10,17)+1
gc=c(1,2,3,8,12,15)}
if(m==30)
{gc=c(0,1,2,5,14,24)+1}

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