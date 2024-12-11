# average variance metrics for one CYCLIC design 
rm(list=ls())
library(MASS)
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
ibloc=c(1,2 ,8, 9,11,   13)}
if(nm==30)
{m=30
nc=6
ibloc=c(1, 2, 3, 6, 15, 25)}

# initial block in coordinates and spacing
gs=gc2s(m,nc,ibloc) 

# generate cyclic kXt design
abmat=gencabd(m,nc,ibloc) 
# number of cyclic sets
tknum(m,nc)
# calculate Ac  
gac=gvabd(m,nc,abmat)
cat("Ac ",round(gac,4))
