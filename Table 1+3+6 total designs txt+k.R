# Table 1+3+6: Square Array Designs 
# generate designs for a t X t square array with k controls
# take m=t and k=nc
rm(list=ls())
library(MASS)
library(tictoc)
library(numbers)
library(gtarray)

tic()

# setting
m=10
nc=4
num=choose(m,nc) # total number of designs
ntot=m^2  # total number of plots
ndes=nc*m  # total number of control plots
ntest=ntot-ndes # total number of test line plots 
ntreat=nc+ntest  # total number of treatments

# basic X and Z matrices
xmat=grxmat(m,nc)
zmat0=grzmat0(m)

# set up coordinates for designs 
gt=t(combn(m,nc))
#matrix of results for coordinates, spacing, and metrics Acc, Act, Att
gtot=matrix(0,num,(2*nc+3)) 
# counters
nprint=1000
ncon=0 

# cycle over all designs
itc=1
for(id in 1:num)
{
# design coordinates and spacing
gcd=gt[id,] #coordinates
gs=gc2s(m,nc,gcd)

# rearrange basic zmat0
zmat=grnewzmat(m,nc,gcd,zmat0)
cmat=grcmat(m,xmat,zmat)
crank=qr(cmat)$rank

# design not connected
if(crank<(ntreat-1))
{ncon=ncon+1
itc=itc+1}

# design connected
if(crank==(ntreat-1))
{
actvec=grmets(m,nc,ntest,ntreat,cmat)
gtot[itc,]=c(gcd,gs,actvec)
if(itc==floor(itc/nprint)*nprint)
{print(itc)}
itc=itc+1
}
}

# order results by Act and hence Att with rounding
gtot1=round(gtot[order(gtot[,(2*nc+2)]), ],6)

# classes for act att 
gactvec=gtot1[,(2*nc+2)]
gattvec=gtot1[,(2*nc+3)]
g1=as.data.frame(table(gactvec))
g2=as.data.frame(table(gattvec))
cbind(g1,g2)
dim(g1)

# number of cyclic sets
tknum(m,nc)

# Act and Att optimal designs
gtot2=subset(gtot1,gtot1[,(2*nc+2)]>0)
glim=min(gtot2[,(2*nc+2)])
gopt=subset(gtot2,gtot2[,(2*nc+2)]<=glim+0.00001)
a=dim(gopt)
#gopt[1:a[1],]
gopt[1:10,]

toc()
