#' Converts coordinates to spacing for cyclic designs
#'
#' @param m    Number of treatments
#' @param nc   Number of rows/controls
#' @param gc   Coordinates -  vector
#'
#' @return gs Spacing
#'
gc2s=function(m,nc,gc)
{if(length(gc) != nc)
{print("invalid nc or length of length(gc)")
  gs=0}else{gs=c(gc[2]-gc[1])
  for(i in 3:nc)
  {gs=c(gs,c(gc[i]-gc[(i-1)]))}
  gs=c(gs,m+gc[1]-gc[nc])}

  return(gs)
}

#' Converts spacing to coordinates for cyclic designs
#'
#' @param m      Number of treatments
#' @param i1     Starting coordinate in first row
#' @param space  Spacing - vector
#' @param nc     Number of rows/controls
#'
#' @return gc    Coordinates
#'
gs2c=function(m,nc,i1,space)
{
  if(sum(space) != m)
  {print("invalid space or m")
    gc=0}else{
      gc0=as.matrix(c(cumsum(space[1:(nc-1)])))
      gc=c(i1,gc0+i1)
      for(i in 2:nc)
      {if(gc[i]>m){gc[i]=gc[i]%%m}}
    }
  return(gc)}

#' Generate a cyclic auxiliary block design from an initial block
#'
#' @param m   Number of treatments
#' @param nc. Number of controls
#' @param ib  Initial block - vector
#'
#' @return abd Augmented block design. Otherwise
#'             0 for an invalid entry
#'
gencabd=function(m,nc,ib)
{
  if(length(ib) != nc)
  {print("invalid nc or length of ib")
    abd=0}
  if(max(ib) > m)
  {print("invalid number of treatments in ib")
    abd=0}
  if(length(ib)== nc & max(ib) <= m)
  {abd=matrix(0,nc,m)
  for(i in 1:nc)
  {
    if(ib[i]==1){abd[i,]=seq(1:m)}else
    {abd[i,]=(c(seq(ib[i],m),seq(1:(ib[i]-1))))}
  }
  }
  return(abd)
}

#' Calculates treatment information matrix C and metric Aabd
#' for an auxiliary block design
#'
#' @param m     Number of treatments
#' @param nc    Number of rows
#' @param abmat Auxiliary block design r x t
#'
#' @return.     Metric Aabd  0 if matrix C is not connected
#'
gvabd=function(m,nc,abmat)
{
  # requisite matrices for calculating the N, R, C matrix
  nr=nrow(abmat)
  ncl=ncol(abmat)
  imat=diag(m)
  jmat=matrix(1,m,m)
  ijmat=imat-(jmat/m)
  # calculate incidence matrix nmat and rmat
  rmat=nc*diag(m)
  nmat=matrix(0,m,m)
  for(r in 1:nr)
  {for(j in 1:ncl)
    nmat[j,abmat[r,j]]=1
  }
  # calculate C mat
  cmat=rmat-(nmat)%*%t(nmat)/nc
  crank=qr(cmat)$rank
  cmatmp=ginv(cmat)
  if(crank<(m-1))
  {gac=0}
  if(crank==(m-1))
  {
    # variances
    gac=2*sum(diag(cmatmp%*%ijmat))/(m-1)
  }
  return(gac)
}

#' Finds number of cyclic sets in a cyclic design
#'
#' @param nt  Number of parameters
#' @param k   Number of row/controls
#'
#' @return    Nsum Number of classes
#'
tknum=function(nt,k)
{dt=gacd(nt)$pos
dk=gacd(k)$pos
dtk=intersect(dt,dk)
Nsum=0
for(d in dtk)
{term=eulersPhi(d)*factorial(nt/d)/(factorial(k/d)*factorial((nt-k)/d))
Nsum=Nsum+term}
Nsum=Nsum/nt
return(Nsum)
}

#' Factors of an integer for calculating
#' the number of cyclic sets
#'
#' @param x Integer
#'
#' @return factors Factors of the integer x
#'
gacd <- function(x) {
  x <- as.integer(x)
  div <- seq_len(abs(x))
  factors <- div[x %% div == 0L]
  factors <- list(neg = -factors, pos = factors)
  return(factors)
}

#' Equivalence classes for t-(k-1) <= 10 ONLY
#'
#' @param nc   Number of rows/controls
#' @param set  Set of designs with the same spacing
#'
#' @return ac  List of equivalence classes (not permuted)
#'
gclass=function(nc,set)
{
  gspace=set[,(nc+1):(2*nc)]
  for(i in 1:nrow(gspace))
  {gspace[i,]=sort(gspace[i,])}
  tvec=c(rev(10^seq(1,(nc-1))),1)
  ac=as.matrix(vec_count(gspace%*%tvec), sort = "key")
  return(ac)
}
