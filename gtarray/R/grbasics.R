#' Calculates X matrix for square array
#'
#' @param m     Dimension of square array m^2
#' @param nc    Number of controls
#'
#' @return xmat X matrix for square array
#'
grxmat=function(m,nc)
{
  # nc controls
  onevec=matrix(1,m,1)
  xmatc1=diag(nc) %x% onevec
  xmatc2=matrix(0,nc*m,m*(m-nc))
  xmatc=cbind(xmatc1,xmatc2)
  # testlines
  xmatt1=matrix(0,m*(m-nc),nc)
  xmatt2=diag(m*(m-nc))
  xmatt=cbind(xmatt1,xmatt2)
  xmat=rbind(xmatc,xmatt)
  return(xmat)
}

#' Calculates Z matrix for square array
#'
#' @param m       Dimension of square array m^2
#'
#' @return zmat   Basic Z matrix 2m x2m
#'
grzmat0=function(m)
{
  # Z ROWS
  onevec=matrix(1,m,1)
  imat=diag(m)
  zmatr=onevec %x% imat
  # Z COLS
  zmatc=imat
  cv=c(m,seq(1:(m-1))) # cycling for diagonals
  gm=imat
  for(i in 1:(m-1))
  {gm=gm[,cv]
  zmatc=rbind(zmatc,gm)
  }
  # Z matrix
  zmat0=cbind(zmatr,zmatc)
  return(zmat0)
}

#' Arrange Zmat0 according to the coordinates of the first row
#'
#' @param m     Dimension of square array m^2
#' @param nc    Number of controls
#' @param gc    Coordinates of controls in first row of the array
#' @param zmat0 Basic Z matrix
#'
#' @return zmat Z matrix for coordinates gc
#'
grnewzmat=function(m,nc,gc,zmat0)
{
  bloc=0
  for(i in gc)
  {a=(i-1)*m+1
  b=i*m
  bloc=as.matrix(c(bloc,seq(a,b)))}
  bloc=bloc[2:(nc*m+1)]
  zmatcon=zmat0[bloc,]
  zmattest=zmat0[-bloc,]
  zmat=rbind(zmatcon,zmattest)
  return(zmat)
}

#' Calculates treatment information matrix C=X^t(I-Z(Z^TZ)^-)X of
#' a square array design
#'
#' @param m      Dimension of square array m^2
#' @param xmat   X matrix for treatments
#' @param zmat   Z matrix for rows and columns
#'
#' @return cmat  C matrix for the square array
#'
grcmat=function(m,xmat,zmat)
{
  intmat=diag(m^2)-zmat%*% ginv(t(zmat) %*% zmat) %*% t(zmat)
  cmat=t(xmat) %*% intmat %*% xmat
  return(cmat)
}

#' Calculates A metrics for a square array design
#'
#' @param m      Dimension of square array m^2
#' @param nc     Number of controls
#' @param cmat   C matrix of square array
#'
#' @return gvct  Variances for pairwise controls - test lines
#'
grmets=function(m,nc,ntest,ntreat,cmat)
{
  cmatmp=ginv(cmat)
  actvec=matrix(0,1,3)
# Acc
  acc=2/m

  #Act
  act=matrix(0,nc*ntest,1)
  it=1
  for(i in 1:nc)
  {for(j in (nc+1):(ntest+nc))
  {cvec=matrix(0,ntreat,1)
  cvec[i]=1
  cvec[j]=-1
  act[it]=t(cvec) %*% cmatmp %*% cvec
  it=it+1
  }
  }
#Att
 cmattt=cmatmp[(nc+1):ntreat,(nc+1):ntreat]
 imat=diag(ntest)
 jmat=matrix(1,ntest,ntest)
 ijmat=imat-(jmat/ntest)
 att=2*sum(diag(cmattt%*%ijmat))/(ntest-1)
#results
    actvec=c(acc,mean(act),mean(att))
    return(actvec)
  }

#' Convert an auxiliary block design to the treatment
#' information matrix of the corresponding  square array design
#'
#' @param m   Number of blocks m
#' @param nc  Number of rows
#' @param abd Auxiliary block design
#'
#' @return   cmat C matrix of corresponding square array design
#'
gab2sq=function(m,nc,abd)
{
  #convert to square array
  sqmat=matrix(0,m,m)
  for(r in 1:nc)
  {for(j in 1:m)
    sqmat[j,abd[r,j]]=r}
  #treatment matrix by ROWS
  tobs=c(t(sqmat))
  xmatc=matrix(0,ntot,nc)
  for(i in 1:ntot)
  {for(j in 1:nc)
  {if(tobs[i]==j)
  {xmatc[i,j]=1}
  }
  }
  # test lines matrix by ROWS
  xmatt=matrix(0,ntot,ntest)
  it=1
  for(i in 1:ntot)
  {if(tobs[i]==0)
  {xmatt[i,it]=1
  it=it+1}
  }
  # full xmatrix
  xmat=cbind(xmatc,xmatt)
  # row col matrix
  onem=matrix(1,m,1)
  zrows=diag(m)%x%onem
  zcols=onem%x%diag(m)
  zmat=cbind(zrows,zcols)
  ztz=t(zmat)%*%zmat
  gztz=ginv(ztz)
  zgztz=zmat%*%gztz%*%t(zmat)
  # cmat
  intmat=diag(m^2)-zmat%*% ginv(t(zmat) %*% zmat) %*% t(zmat)
  cmat=t(xmat) %*% intmat %*% xmat
  return(cmat)
}

