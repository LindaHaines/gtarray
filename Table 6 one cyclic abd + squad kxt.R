# average variance metrics for a cyclic auxiliary 
# block design and corresponding square array design
rm(list=ls())
library(MASS)
library(numbers)
library(gtarray)

# examples
nm=16
# setting
if(nm==12)
{m=12
nc=3
ibloc=c(1,4,8)}
if(nm==16)
{m=16
nc=6
ibloc=c(1,2,8,9,11,13)}
if(nm==30)
{m=30
nc=6
ibloc=c(1,2,3,6,15,25)}

# Auxiliary Block Design
 
# generate cyclic kXt design
gs=gc2s(m,nc,ibloc) 
abmat=gencabd(m,nc,ibloc)
# number of cyclic sets
tknum(m,nc)
# calculate Ac  
gac=gvabd(m,nc,abmat)
 
if(gac==0) {print("not connected")
  print(gac)
  print(ibloc)
  print(gs)}

# Square Array Design 

# numbers of treatments and plots
num=choose(m,nc) 
ntot=m^2
ndes=nc*m
ntest=ntot-ndes
ntreat=nc+ntest 

# generate square array design from auxiliary block design
if(gac > 0)
{ print("connected")
# convert abd to square array  Acc Act Att
cmat=gab2sq(m,nc,abmat)
# generate A metrics
actvec=grmets(m,nc,ntest,ntreat,cmat)
}

# setting and Ac, Acc Act Att
cat("t=",m,"k=",nc,"initial block", ibloc,"\n")
cat("Ac", gac,"\n")
cat("Acc Act Att", actvec)
