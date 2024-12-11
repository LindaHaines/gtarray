# columns for Table 5 k=4
rm(list=ls())
library(MASS)
library(gtarray)

# setting
nc=4
for(m in 14:26)
{if(m==14)
{gc=c(0,1,4,6)+1}
if(m==15 | m==17|m==18)
{gc=c(0,1,3,7)+1}
if(m==16)
{gc=c(0,1,3,12)+1}
if(m==19)
{gc=c(0,1,3,8)+1}
if(m==20)
{gc=c(0,1,3,14)+1}
if(m==21 | m==22)
{gc=c(0,1,3,9)+1}
if(m==23 | m==24|m==25)
{gc=c(0,1,3,10)+1}
if(m==26)
{gc=c(0,1,3,11)+1}
  
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