# calculate Ac, Acc, Act, Att for a BIBD
rm(list=ls())

# setting
m=16
nc=6
lam=nc*(nc-1)/(m-1)
m1=m*(m-nc) # number of plots for test lines

# formulae
ac=2*nc/(m*lam)
acc=2/m;
att=2+4*(m-1)*(m-nc)/((m1-1)*(nc-1))
act=1+(2*nc+lam)/(m*lam)
actvec=c(acc,act,att)

# results 
# setting and Ac, Acc Act Att
cat("t=",m,"k=",nc,"lambda", lam,"\n")
cat("Ac", round(ac,4),"\n")
cat("Acc Act Att", round(actvec,4))

