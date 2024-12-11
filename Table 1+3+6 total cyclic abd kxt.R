# Table 1+3+6: Auxiliary Block Designs
# generate all designs for all cyclic k X t auxiliary block designs
# for ease of programming take m=t and k=nc
# NOTE: package ibd could also be used here
rm(list=ls())
library(MASS)
library(tictoc)
library(numbers) 
library(vctrs) 
library(gtarray)

tic()
# setting
m=9
nc=3
num=choose(m,nc) # total number of designs

# coordinates for cycling through all possible designs 
gt=t(combn(m,nc))
ngt=nrow(gt)

# counter for number of cycles
nprint=10000

# set up matrix of results for the designs
# coordinates,spacings and metric Ac
gtot=matrix(0,num,(2*nc+1))

# cycle over all designs
itc=1
for(ndes in 1:ngt)
{
# generate cyclic auxiliary block design from initial block
  ib=gt[ndes,]
  abmat=gencabd(m,nc,ib)  
# calculate Ac and save 
  gac=gvabd(m,nc,abmat)
  gs=gc2s(m,nc,gt[ndes,])
  gtot[ndes,]=c(gt[ndes,],gs,gac)
# monitor cycles
  if(itc==floor(itc/nprint)*nprint)
  {print(itc)}
  itc=itc+1
}

# total number of cyclic sets
tknum(m,nc)

# Ac and number of equivalence classes
gvec=round(gtot,6)
g1=as.data.frame(table(gvec[,(2*nc+1)]))
g1

# Ac-optimal designs: coordinates, spacing, Ac metric
gtot1=round(gtot[order(gtot[,(2*nc+1)]), ],10)
gtot2=subset(gtot1,gtot1[,(2*nc+1)]>0)
glim=min(gtot2[,(2*nc+1)])
gopt=subset(gtot2,gtot2[,(2*nc+1)]<=glim+0.00001)
dim(gopt)
gopt[1:10,]

# CAUTION! 
# Equivalence classes with no permuting 
# Only for k=3 and t-(k-1) < = 10 with t-(k-1) = 10 indicated as "max"
if(nc ==3 & m-(nc-1) <= 10 )
{dc=cumsum(as.matrix(g1[,2]))
nd=length(dc)
gequiv=gtot1[1:dc[1],]
print(gclass(nc,gequiv))
for(i in 1:(nd-1))
{
  gequiv=gtot1[(dc[i]+1):dc[i+1],]
  print(gclass(nc,gequiv))
  mx=max(gequiv[,(nc+1):(2*nc)])
  if(mx >9){print(paste0("max  ",max(gequiv[,(nc+1):(2*nc)])))}
}}

toc()
