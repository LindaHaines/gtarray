# columns for Table 5 k=5
rm(list=ls())
library(MASS)
library(gtarray)

# setting
nc=5
for(m in 17:30)
{if(m==17)
{gc=c(0,1,3,7,8)+1}
if(m==18)
{gc=c(0,1,2,5,11)+1}
if(m==19)
{gc=c(0,1,2,6,9)+1}
if(m==20)
{gc=c(1,2,4,8,13)}
#gc=c(0,1,2,5,14)+1}
if(m==21)
{gc=c(1,2,5,15,17)}
if(m==22)
{gc=c(0,1,3,7,12)+1}
if(m==23)
{gc=c(0,1,3,8,14)+1}
if(m==24)
{gc=c(0,1,3,11,20)+1}
if(m==25)
{gc=c(0,1,3,15,21)+1}
if(m==26)
{gc=c(0,1,3,7,12)+1}  
if(m==27)
{gc=c(0,1,3,7,18)+1}
if(m==28)
{gc=c(0,1,3,13,24)+1}
if(m==29)
{gc=c(0,1,3,7,21)+1}
if(m==30)
{gc=c(0,1,3,12,25)+1} 
  
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