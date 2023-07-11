################################################
################## R PACKAGES ##################
################################################

library(MASS)
library(fda)
library(gmp)
library(kernlab)
library(SwarmSVM)
library(factoextra)
library(dplyr)
library(cluster)
library(xtable)

###############################################
################## FUNCTIONS ##################
###############################################

## All the functions created for this paper are defined in this R script.

######################
### Epigraph Index ###
######################

## The epigraph index of a curve x (EI) is one minus the proportion of curves in the sample that lie above x.
# This function computes the EI of a bunch of B curves.

EI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.min <- curves[i,]
    for (j in 1:B){
      inpoints <- sum((curves[j,] >= env.min))
      if (inpoints == lengthcurves)
      {index[i] <- index[i]+1}
    }
  }
  index <- index/B
  return (1-index)
}

#######################
### Hypograph Index ###
#######################

## The hypograph index of a curve x (HI) is the proportion of curves in the sample that lie below x.
# This function computes the HI of a bunch of B curves.

HI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.max <- curves[i,]
    for (j in 1:B){
      inpoints <- sum(curves[j,]<=env.max)
      if (inpoints == lengthcurves)
      {index[i] <- index[i]+1}
    }
  }
  index <- index/B
  return (index)
}

##################################
### Generalized Epigraph Index ###
##################################

## The generalized epigraph index of a curve x (MEI) is 
## one minus the "proportion of time" that the curves of the sample are above x.
# This function computes the MEI of a bunch of B curves.

MEI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.min <- curves[i,]
    for (j in 1:B){
      inpoints <- sum(curves[j,] >= env.min)
      index[i] <- index[i]+inpoints
    }
  }
  index <- index/(B*lengthcurves)
  return (1-index)
}

###################################
### Generalized Hypograph Index ###
###################################

## The generalized hypograph index of a curve x (MHI) is 
## the "proportion of time" that the curves of the sample are below x.
# This function computes the MHI of a bunch of B curves.

MHI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.max <- curves[i,]
    for (j in 1:B){
      inpoints <- sum(curves[j,]<=env.max)
      index[i] <- index[i]+inpoints
    }
  }
  index <- index/(B*lengthcurves)
  return (index)
}

######################
### Multivariate Epigraph Index ###
######################

## The epigraph index of a multivariate curve x (mulEI) is one minus the 
# proportion of curves in the sample that lie above x (considering x as a multivariate curve).
# This function computes the mulEI of a bunch of B curves.

mulEI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  D <- dim(curves)[3]
  for (i in 1:B){
    index[i] <- 0
    inpoints_k <- array(rep(NaN,B*lengthcurves*D),dim=c(B,lengthcurves,D))
    for (k in 1:D){
      env.min <- curves[i,,k]
      inpoints <- c()
      for (j in 1:B){
        inpoints<- rbind(inpoints,curves[j,,k] >= env.min)
      }
      inpoints_k[,,k]<-inpoints
    }
    inpoints_mult<-inpoints_k[,,1];
    for (l in 2 : D){
      inpoints_mult<- inpoints_mult*inpoints_k[,,l]
    }
    for (m in 1:B){
      inpoints_f <- sum(inpoints_mult[m,])
      if (inpoints_f == lengthcurves)
      {index[i] <- index[i]+1}
    }
  }
  index <- index/B
  return (1-index)
}

## The hipograph index of a multivariate curve x (mulHI) is the 
# proportion of curves in the sample that lie below x (considering x as a multivariate curve).
# This function computes the mulHI of a bunch of B curves.

mulHI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  D <- dim(curves)[3]
  for (i in 1:B){
    index[i] <- 0
    inpoints_k <- array(rep(NaN,B*lengthcurves*D),dim=c(B,lengthcurves,D))
    for (k in 1:D){
      env.min <- curves[i,,k]
      inpoints <- c()
      for (j in 1:B){
        inpoints<- rbind(inpoints,curves[j,,k] <= env.min)
      }
      inpoints_k[,,k]<-inpoints
    }
    inpoints_mult<-inpoints_k[,,1];
    for (l in 2 : D){
      inpoints_mult<- inpoints_mult*inpoints_k[,,l]
    }
    for (m in 1:B){
      inpoints_f <- sum(inpoints_mult[m,])
      if (inpoints_f == lengthcurves)
      {index[i] <- index[i]+1}
    }
  }
  index <- index/B
  return (index)
}

## The generalized hipograph index of a multivariate curve x (mulMHI) is one minus the 
# "proportion of time" that curves in the sample lie below x (considering x as a multivariate curve).
# This function computes the mulMHI of a bunch of B curves.

mulMHI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  D <- dim(curves)[3]
  for (i in 1:B){
    index[i] <- 0
    inpoints_k <- array(rep(NaN,B*lengthcurves*D),dim=c(B,lengthcurves,D))
    for (k in 1:D){
      env.min <- curves[i,,k]
      inpoints <- c()
      for (j in 1:B){
        inpoints<- rbind(inpoints,curves[j,,k] <= env.min)
      }
      inpoints_k[,,k]<-inpoints
    }
    inpoints_mult<-inpoints_k[,,1];
    for (l in 2 : D){
      inpoints_mult<- inpoints_mult*inpoints_k[,,l]
    }
    for (m in 1:B){
      inpoints_f <- sum(inpoints_mult[m,])
      index[i] <- index[i]+inpoints_f
    }
  }
  index <- index/(B*lengthcurves)
  return (index)
}

## The generalized epigraph index of a multivariate curve x (mulMEI) is one minus the 
# "proportion of time" that curves in the sample lie above x (considering x as a multivariate curve).
# This function computes the mulMEI of a bunch of B curves.

mulMEI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  D <- dim(curves)[3]
  for (i in 1:B){
    index[i] <- 0
    inpoints_k <- array(rep(NaN,B*lengthcurves*D),dim=c(B,lengthcurves,D))
    for (k in 1:D){
      env.min <- curves[i,,k]
      inpoints <- c()
      for (j in 1:B){
        inpoints<- rbind(inpoints,curves[j,,k] >= env.min)
      }
      inpoints_k[,,k]<-inpoints
    }
    inpoints_mult<-inpoints_k[,,1];
    for (l in 2 : D){
      inpoints_mult<- inpoints_mult*inpoints_k[,,l]
    }
    for (m in 1:B){
      inpoints_f <- sum(inpoints_mult[m,])
      index[i] <- index[i]+inpoints_f
    }
  }
  index <- index/(B*lengthcurves)
  return (1-index)
}

##########################
### funspline function ###
##########################

## This function fits a smooth B-spline basis for a bunch of curves X and
## calculates first and second derivatives for each curve.

funspline <- function(X,t,basis){
  
  # INPUT
  # X <- functional data object
  # t <- sequence for creating the basis
  # basis <- create.bspline.basis() object
  
  ys = smooth.basis(argvals=t, y=t(X), fdParobj=basis) #data obtained when applying a smoothed B-spline basis
  smooth <- t(eval.fd(t,ys$fd,0)) #smoothed data
  deriv <- t(eval.fd(t,ys$fd,1)) #first derivatives
  deriv2 <- t(eval.fd(t,ys$fd,2)) #second derivatives
  
  res <- list("spline"=ys, "smooth"=smooth, "deriv"=deriv, "deriv2"=deriv2) #object containing the data and 1st and 2nd derivatives
  return(res)
}

####################
### ind function ###
####################

## This function generates a dataframe containing EI, HI, MEI and MHI
## for each curve on the original data, first and second derivatives.

ind <- function(X,t,rangeval=c(0,1),nbasis=40,norder=4){
  
  # INPUT
  # X <- functional data object containing the different groups that can be originally observed
  #      (1 to x rows correspond to the observations for the first group, x+1 to y correspond to the second group and so on)
  # t <- sequence for creating the basis
  # rangeval <- a numeric vector of length 2 defining the interval over which the
  #             functional data object can be evaluated. By default, we consider the interval [0,1]
  # nbasis <- number of basis functions 
  # norder <- order of b-splines
  
  basisobj <- create.bspline.basis(rangeval=rangeval,norder=norder,nbasis=nbasis) #b-spline basis (order 3 polynomials when norder = 4)
  
  # dta <- X #original data
  
  if(length(dim(X))==3){
    # Multivariate functional data
    
    N <- dim(X)[1]
    P <- dim(X)[2]
    K <- dim(X)[3]
    
    dta <- array(rep(NaN,N*P),dim=c(N,P,K))
    ddta <- array(rep(NaN,N*P),dim=c(N,P,K))
    d2dta <- array(rep(NaN,N*P),dim=c(N,P,K))
    for(k in 1:K){
      dta_funspline <- funspline(X[,,k],t,basisobj)
      dta[,,k] <- dta_funspline$smooth #smoothed data
      ddta[,,k] <- dta_funspline$deriv #first derivative
      d2dta[,,k] <- dta_funspline$deriv2 #second derivative
    }  
    # Applying the indexes to the origial data
    dtaEI <- mulEI(dta) 
    dtaHI <- mulHI(dta)
    dtaMEI <- mulMEI(dta)
    dtaMHI <- mulMHI(dta)
    
    # Applying the indexes to the data first derivatives
    ddtaEI <- mulEI(ddta)
    ddtaHI <- mulHI(ddta)
    ddtaMEI <- mulMEI(ddta)
    ddtaMHI <- mulMHI(ddta)
    
    # Applying the indexes to the data second derivatives
    d2dtaEI <- mulEI(d2dta)
    d2dtaHI <- mulHI(d2dta)
    d2dtaMEI <- mulMEI(d2dta)
    d2dtaMHI <- mulMHI(d2dta)
    
    # New multivariate data set
    ind.data <- data.frame(dtaEI,dtaHI,dtaMEI,dtaMHI,ddtaEI,ddtaHI,ddtaMEI, 
                           ddtaMHI,d2dtaEI,d2dtaHI,d2dtaMEI,d2dtaMHI)
  } else if(length(dim(X))==2){
    # Univariate functional data
    
    dtaX <- funspline(X,t,basisobj)
    dta <- dtaX$smooth # smoothed data coefs
    ddta <- dtaX$deriv #first derivatives
    d2dta <- dtaX$deriv2 #second derivatives
    
    dtaEI <- EI(dta) 
    dtaHI <- HI(dta)
    dtaMEI <- MEI(dta)
    dtaMHI <- MHI(dta)
    
    # Applying the indexes to the data first derivatives
    ddtaEI <- EI(ddta)
    ddtaHI <- HI(ddta)
    ddtaMEI <- MEI(ddta)
    ddtaMHI <- MHI(ddta)
    
    # Applying the indexes to the data second derivatives
    d2dtaEI <- EI(d2dta)
    d2dtaHI <- HI(d2dta)
    d2dtaMEI <- MEI(d2dta)
    d2dtaMHI <- MHI(d2dta)
    
    # New multivariate data set
    ind.data <- data.frame(dtaEI,dtaHI,dtaMEI,dtaMHI,ddtaEI,ddtaHI,ddtaMEI, 
                           ddtaMHI,d2dtaEI,d2dtaHI,d2dtaMEI,d2dtaMHI)
  } else{
    print("Non valid data dimension")
    ind.data <- "Non valid data dimension"
  }
  
  
  return(ind.data)
}

################################################
### ind function (only generalized indexes) ###
###############################################

## This function generates a dataframe containing EI, HI, MEI and MHI
## for each curve on the original data, first and second derivatives.

indM <- function(X,t,rangeval=c(0,1),nbasis=40,norder=4){
  
  # INPUT
  # X <- functional data object containing the different groups that can be originally observed
  #      (1 to x rows correspond to the observations for the first group, x+1 to y correspond to the second group and so on)
  # t <- sequence for creating the basis
  # rangeval <- a numeric vector of length 2 defining the interval over which the
  #             functional data object can be evaluated. By default, we consider the interval [0,1]
  # nbasis <- number of basis functions 
  # norder <- order of b-splines
  
  basisobj <- create.bspline.basis(rangeval=rangeval,norder=norder,nbasis=nbasis) #b-spline basis (order 3 polynomials when norder = 4)
  
  # dta <- X #original data
  
  if(length(dim(X))==3){
    # Multivariate functional data
    
    N <- dim(X)[1]
    P <- dim(X)[2]
    K <- dim(X)[3]
    
    dta <- array(rep(NaN,N*P),dim=c(N,P,K))
    ddta <- array(rep(NaN,N*P),dim=c(N,P,K))
    d2dta <- array(rep(NaN,N*P),dim=c(N,P,K))
    for(k in 1:K){
      dta_funspline <- funspline(X[,,k],t,basisobj)
      dta[,,k] <- dta_funspline$smooth #smoothed data
      ddta[,,k] <- dta_funspline$deriv #first derivative
      d2dta[,,k] <- dta_funspline$deriv2 #second derivative
    }  
    # Applying the indexes to the origial data
    # dtaEI <- mulEI(dta) 
    # dtaHI <- mulHI(dta)
    dtaMEI <- mulMEI(dta)
    dtaMHI <- mulMHI(dta)
    
    # Applying the indexes to the data first derivatives
    # ddtaEI <- mulEI(ddta)
    # ddtaHI <- mulHI(ddta)
    ddtaMEI <- mulMEI(ddta)
    ddtaMHI <- mulMHI(ddta)
    
    # Applying the indexes to the data second derivatives
    # d2dtaEI <- mulEI(d2dta)
    # d2dtaHI <- mulHI(d2dta)
    d2dtaMEI <- mulMEI(d2dta)
    d2dtaMHI <- mulMHI(d2dta)
    
    # New multivariate data set
    ind.data <- data.frame(dtaMEI,dtaMHI,ddtaMEI,ddtaMHI,d2dtaMEI,d2dtaMHI)
  } else if(length(dim(X))==2){
    # Univariate functional data
    
    dtaX <- funspline(X,t,basisobj)
    dta <- dtaX$smooth # smoothed data coefs
    ddta <- dtaX$deriv #first derivatives
    d2dta <- dtaX$deriv2 #second derivatives
    
    dtaEI <- EI(dta) 
    dtaHI <- HI(dta)
    dtaMEI <- MEI(dta)
    dtaMHI <- MHI(dta)
    
    # Applying the indexes to the data first derivatives
    ddtaEI <- EI(ddta)
    ddtaHI <- HI(ddta)
    ddtaMEI <- MEI(ddta)
    ddtaMHI <- MHI(ddta)
    
    # Applying the indexes to the data second derivatives
    d2dtaEI <- EI(d2dta)
    d2dtaHI <- HI(d2dta)
    d2dtaMEI <- MEI(d2dta)
    d2dtaMHI <- MHI(d2dta)
    
    # New multivariate data set
    ind.data <- data.frame(dtaEI,dtaHI,dtaMEI,dtaMHI,ddtaEI,ddtaHI,ddtaMEI, 
                           ddtaMHI,d2dtaEI,d2dtaHI,d2dtaMEI,d2dtaMHI)
  } else{
    print("Non valid data dimension")
    ind.data <- "Non valid data dimension"
  }
  
  
  return(ind.data)
}

######################
### valid function ###
######################

## This function applies the clustering evaluation criteria considered during
## the memory (Purity, F-measure, Rand Index(RI)).

valid <- function(true_labels, clusters){
  
  # INPUT
  # true_labels <- vector containing the correct classification
  # clusters <- vector containing the classification obtained by a clustering method
  
  if (is.integer(true_labels)) 
    true_labels = as.numeric(true_labels)
  if (is.integer(clusters)) 
    clusters = as.numeric(clusters)
  if (!is.vector(true_labels) || !is.numeric(true_labels)) 
    stop("true_labels should be a numeric vector")
  if (!is.vector(clusters) || !is.numeric(clusters)) 
    stop("clusters should be a numeric vector")
  if (length(true_labels) != length(clusters)) 
    stop("The length of the true_labels vector should be equal to the length of the clusters vector")
  
  tbl = table(clusters, true_labels) #contingency table
  conv_df = as.data.frame.matrix(tbl)
  
  # Purity
  tmp_pur = apply(conv_df, 1, max)
  res_purity = sum(tmp_pur)/length(true_labels)
  
  # True positives(tp), false positives(fp), true negatives(tn) and false negatives(fn)
  # (needed for calculating RI)
  tp_plus_fp = sum(gmp::asNumeric(gmp::chooseZ(rowSums(conv_df), 2)))
  tp_plus_fn = sum(gmp::asNumeric(gmp::chooseZ(colSums(conv_df), 2)))
  tp = sum(gmp::asNumeric(gmp::chooseZ(as.vector(as.matrix(conv_df)), 2)))
  fp = tp_plus_fp - tp
  fn = tp_plus_fn - tp
  tn = gmp::asNumeric(gmp::chooseZ(sum(as.vector(as.matrix(conv_df))), 2)) - tp - fp - fn
  
  # Precision and recall (needed for calculating Fmeasure)
  prec = tp/(tp + fp) # Precision
  rec = tp/(tp + fn) # Recall
  
  # (Purity, Fmeasure, RI)
  res <- c(round(res_purity, 4), round(2 *((prec * rec)/(prec + rec)), 4),
           round((tp + tn)/(tp + fp + fn + tn), 4))
  names(res) <- c("Purity", "Fmeasure", "RI")
  
  return(as.table(res))
}

#############################
### kmeans_mahal function ###
#############################

## This function computes k-means multivariate clustering with Mahalanobis distance.

kmeans_mahal <- function(dta,clus){
  
  #INPUT
  # dta <- dataframe containing data for applying k-means clustering
  # clus <- number of clusters
  
  x <- as.matrix(dta)
  c <- chol(var(x))
  y <- x %*% solve(c)
  km <- kmeans(dta,centers=clus,iter.max=1000,nstart=100)$cluster #vector containing the clustering partition
  return(km)
  
}

#############################
### clustInd_hierarch_aux ###
#############################

clustInd_hierarch_aux <- function(dtaset, VARS, clust_method, dist, clus, true_labels){
  
  t0 <- Sys.time()
  d <-dist(dtaset[,VARS], method = dist)
  met <- hclust(d, method = clust_method)
  clus_part <- cutree(met, k = clus)
  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))
  if(missing(true_labels)){
    res <- list("cluster"=clus_part, "time"=as.numeric(t))
  } else{
    valid <- valid(true_labels, clus_part)
    res <- list("cluster"=clus_part, "valid"=valid, "time"=as.numeric(t))
  }
  return(res)
}

#########################
### clustInd_hierarch ###
#########################

## This function applies hierarchical clustering methods to a data set obtained from applying indexes to
## data composed by a number equal to clus bunches of curves and also applies validation techniques to the 
## obtained clustering results. Execution times are also computed.

clustInd_hierarch <- function(X_ind,t,rangeval=c(0,1),nbasis=40,norder=4,dist="euclidean",clus=2,true_labels){
  
  #INPUT
  # X <- functional data object containing the different groups that can be originally observed
  # X_ind <- ind(X,t,rangeval,nbasis,norder) #multivariate dataframe (applying indexes)
  
  #      (1 to x rows correspond to the observations for the first group, x+1 to y correspond to the second group and so on)
  # t <- sequence for creating the basis
  # rangeval <- a numeric vector of length 2 defining the interval over which the functional data object can be evaluated
  # nbasis <- number of basis functions
  # norder <- order of b-splines
  # dist <- distance for hierarchical clustering ('euclidean' or 'manhattan')
  # clus <- number of clusters
  # true_labels <- vector containing the correct classification
  
  dtaset <- X_ind
  
  distances <- c("euclidean","manhattan")
  
  if(is.na(pmatch(dist,distances)))
    stop("invalid distance method")
  
  
  # generalized indexes
  VARS1 <- c("dtaMEI","dtaMHI")
  VARS2 <- c("ddtaMEI","ddtaMHI")
  VARS3 <- c("d2dtaMEI","d2dtaMHI")
  VARS4 <- c(VARS1,VARS2)
  VARS5 <- c(VARS1,VARS3)
  VARS6 <- c(VARS2,VARS3)
  VARS7 <- c(VARS4,VARS3)
  
  clust_meth <- c("single","complete","average","centroid","ward.D2")
  VARS <- c("VARS1","VARS2","VARS3","VARS4","VARS5","VARS6","VARS7")
  names <- c("._.MEIMHI",".d.MEIMHI",".d2.MEIMHI","._d.MEIMHI","._d2.MEIMHI",
             ".dd2.MEIMHI","._dd2.MEIMHI")
  
  # APPLYING DIFERENT HIERARCHICAL METHODS
  
  val_indices <- data.frame()
  clusters <- list()
  times <- data.frame()
  for(m in 1:length(clust_meth)){
    names_method <- c()
    name <- c()
    val_indices_i <- c()
    clusters_i <- list()
    times_i <- c()
    for(i in 1:18){
      #if(abs(det(var(dtaset[,eval(parse(text = VARS[i]))])))>1e-5){
        tryCatch({
          c <- clustInd_hierarch_aux(dtaset, eval(parse(text = VARS[i])), clust_meth[m], dist, clus, true_labels)
          name <- paste(clust_meth[m],names[i],"-",dist,sep="")
          names_method <- c(names_method, name)
          clusters_i <- c(clusters_i, list(c$cluster))
          times_i <- rbind(times_i, c$time)
          if(! (missing(true_labels))){
            val_indices_i <- rbind(val_indices_i,c$valid)
          }
        }, error=function(e){}) 
      #}
    }
    
    names(clusters_i) <- names_method
    row.names(times_i) <- names_method
    colnames(times_i)[1] <- "time"
    
    if(! (missing(true_labels))){
      row.names(val_indices_i) <- names_method
      val_indices <- rbind(val_indices, val_indices_i)
    }
    clusters <- c(clusters, clusters_i)
    times <- rbind(times, times_i)
  }
  if(missing(true_labels)){
    res <- list("clusters" = clusters, "time" = times) 
  } else{
    res <- list("val_indices" = as.data.frame(val_indices), "clusters" = clusters, "time" = times) 
  }
  return(res)
}


###########################
### clustInd_kmeans_aux ###
###########################

clustInd_kmeans_aux <- function(dtaset, VARS, dist, clus, true_labels){
  
  t0 <- Sys.time()
  if(dist=="euclidean"){
    clus_part <- kmeans(dtaset[,VARS],centers=clus,iter.max=1000,nstart=100)$cluster
  }else{#mahalanobis
    clus_part <- kmeans_mahal(dtaset[,VARS],clus)
  }
  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))
  if(missing(true_labels)){
    res <- list("cluster"=clus_part, "time"=as.numeric(t))
  } else{
    valid <- valid(true_labels,clus_part)
    res <- list("cluster"=clus_part, "valid"=valid, "time"=as.numeric(t))
  }
  return(res)
}  

#######################
### clustInd_kmeans ###
#######################

## This function applies partitional clustering technique (k-means) to a data set obtained from applying indexes to
## data composed by two bunches of curves and also applies validation techniques to the obtained clustering results.
## Here, execution times are also computed.

clustInd_kmeans <- function(X_ind,t,rangeval=c(0,1),nbasis=40,norder=4,dist='euclidean',clus=2,true_labels){
  
  #INPUT
  # X <- functional data object containing the different groups that can be originally observed
  #      (1 to x rows correspond to the observations for the first group, x+1 to y correspond to the second group and so on)
  # t <- sequence for creating the basis
  # rangeval <- a numeric vector of length 2 defining the interval over which the functional data object can be evaluated
  # nbasis <- number of basis functions
  # norder <- order of b-splines
  # dist <- distance for k-means ('euclidean' or 'mahalanobis') 
  # clus <- number of clusters
  # true_labels <- data frame containing the correct classification
  
  distances <- c("euclidean","mahalanobis")
  if(is.na(pmatch(dist,distances)))
    stop("invalid distance method")
  
  # dtaset <- ind(X,t,rangeval,nbasis,norder) #multivariate data (applying indexes)
  dtaset <- X_ind
  
  # generalized indexes
  VARS1 <- c("dtaMEI","dtaMHI")
  VARS2 <- c("ddtaMEI","ddtaMHI")
  VARS3 <- c("d2dtaMEI","d2dtaMHI")
  VARS4 <- c(VARS1,VARS2)
  VARS5 <- c(VARS1,VARS3)
  VARS6 <- c(VARS2,VARS3)
  VARS7 <- c(VARS4,VARS3)
  
  VARS <- c("VARS1","VARS2","VARS3","VARS4","VARS5","VARS6","VARS7")
  names <- c("._.MEIMHI",".d.MEIMHI",".d2.MEIMHI","._d.MEIMHI","._d2.MEIMHI",
             ".dd2.MEIMHI","._dd2.MEIMHI")
  
  val_indices <- c()
  clusters <- list()
  times <- data.frame()
  names_method <- c()
  name <- c()
  
  for(i in 1:18){
    #if(abs(det(var(dtaset[,eval(parse(text = VARS[i]))])))>1e-5){
      tryCatch({
        c <- clustInd_kmeans_aux(dtaset, eval(parse(text = VARS[i])), dist, clus, true_labels)
        name <- paste("kmeans",names[i],"-",dist,sep="")
        names_method <- c(names_method, name)
        clusters <- c(clusters, list(c$cluster))
        times <- rbind(times, c$time)
        if(!(missing(true_labels))){
          val_indices <- rbind(val_indices,c$valid)
        }
      }, error=function(e){}) 
    #}
  }
  names(clusters) <- names_method
  row.names(times) <- names_method
  colnames(times)[1] <- "time"
  
  if(missing(true_labels)){
    res <- list("clusters" = clusters, "time" = times) 
  }else{
    row.names(val_indices) <- names_method
    res <- list("val_indices" = as.data.frame(val_indices), "clusters" = clusters, "time" = times) 
  }
  return(res)
}

#############################
### clustInd_kkmeans_aux ####
#############################

clustInd_kkmeans_aux <- function(dtaset, VARS, kernel, clus, true_labels){
  
  t0 <- Sys.time()
  met <- kernlab::kkmeans(as.matrix(dtaset[,VARS]), clus, na.action = na.omit, kernel=kernel)
  clus_part <- met@.Data
  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))
  if(missing(true_labels)){
    res <- list("cluster"=clus_part, "time"=as.numeric(t))
  } else{
    valid <- valid(true_labels,clus_part)
    res <- list("cluster"=clus_part, "valid"=valid, "time"=as.numeric(t))
  }
  
  return(res)
} 

##########################
### clustInd_kkmeans ####
#########################

## This function applies kernel k-means to a functional data set composed by two bunches 
## of curves and also applies validation techniques to the obtained clustering results.
## Here, execution times are also computed.

clustInd_kkmeans <- function(X_ind,t,rangeval=c(0,1),nbasis=40,kernel ="rbfdot",norder=4,clus=2,true_labels){
  
  #INPUT
  # X <- functional data object containing the different groups that can be originally observed
  #      (1 to x rows correspond to the observations for the first group, x+1 to y correspond to the second group and so on)
  # t <- sequence for creating the basis
  # rangeval <- a numeric vector of length 2 defining the interval over which the functional data object can be evaluated
  # nbasis <- number of basis functions
  # norder <- order of b-splines
  # kernel <- rbfdot (gaussian), polydot (polynomial) 
  # clus <- number of clusters
  # true_labels <- data frame containing the correct classification
  
  #dtaset <- ind(X,t,rangeval,nbasis,norder) #multivariate data (applying indexes)
  dtaset <- X_ind
  
  # generalized indexes
  VARS1 <- c("dtaMEI","dtaMHI")
  VARS2 <- c("ddtaMEI","ddtaMHI")
  VARS3 <- c("d2dtaMEI","d2dtaMHI")
  VARS4 <- c(VARS1,VARS2)
  VARS5 <- c(VARS1,VARS3)
  VARS6 <- c(VARS2,VARS3)
  VARS7 <- c(VARS4,VARS3)
  
  VARS <- c("VARS1","VARS2","VARS3","VARS4","VARS5","VARS6","VARS7")
  names <- c("._.MEIMHI",".d.MEIMHI",".d2.MEIMHI","._d.MEIMHI","._d2.MEIMHI",
             ".dd2.MEIMHI","._dd2.MEIMHI")  
  val_indices <- c()
  clusters <- list()
  times <- data.frame()
  names_method <- c()
  name <- c()
  
  for(i in 1:18){
    #if(abs(det(var(dtaset[,eval(parse(text = VARS[i]))])))>1e-5){
      tryCatch({
        c <- clustInd_kkmeans_aux(dtaset, eval(parse(text = VARS[i])), kernel, clus, true_labels)
        name <- paste("kkmeans",names[i],"-",kernel,sep="")
        names_method <- c(names_method, name)
        clusters <- c(clusters, list(c$cluster))
        times <- rbind(times, c$time)
        if(!(missing(true_labels))){
          val_indices <- rbind(val_indices,c$valid)
        }
      }, error=function(e){})  
    #}
  }
  names(clusters) <- names_method
  row.names(times) <- names_method
  colnames(times)[1] <- "time"
  
  if(missing(true_labels)){
    res <- list("clusters" = clusters, "time" = times) 
  } else{
    row.names(val_indices) <- names_method
    res <- list("val_indices" = as.data.frame(val_indices), "clusters" = clusters, "time" = times) 
  }
  return(res)
}


#########################
### clustInd_svc_aux ####
#########################

clustInd_svc_aux <- function(dtaset, VARS, cluster.method, clus, true_labels){
  
  t0 <- Sys.time()
  y <- c()
  for(i in 1:(clus-1)){
    y <- c(y, rep(i,as.integer((dim(dtaset)[1]/clus))))
  }
  le <- length(y)
  y <- c(y,rep(clus, dim(dtaset)[1]-le))
  met <- clusterSVM(dtaset[,VARS], y, clus, cluster.method=cluster.method)
  clus_part <- met$label
  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))
  if(missing(true_labels)){
    res <- list("cluster"=clus_part, "time"=as.numeric(t))
  }else{
    valid <- valid(true_labels,clus_part)
    res <- list("cluster"=clus_part, "valid"=valid, "time"=as.numeric(t))
  }
  return(res)
} 

#####################
### clustInd_svc ####
#####################

## This function applies support vector clustering to a functional data set composed by two bunches 
## of curves and also applies validation techniques to the obtained clustering results.
## Here, execution times are also computed.

clustInd_svc <- function(X_ind,t,rangeval=c(0,1),nbasis=40,cluster.method="kmeans",norder=4,clus=2,true_labels){
  
  #INPUT
  # X <- functional data object containing the different groups that can be originally observed
  #      (1 to x rows correspond to the observations for the first group, x+1 to y correspond to the second group and so on)
  # t <- sequence for creating the basis
  # rangeval <- a numeric vector of length 2 defining the interval over which the functional data object can be evaluated
  # nbasis <- number of basis functions
  # norder <- order of b-splines
  # cluster.method <- kmeans or kernkmeans 
  # clus <- number of clusters
  # true_labels <- data frame containing the correct classification
  
  # dtaset <- ind(X,t,rangeval,nbasis,norder) #multivariate data (applying indexes)
  dtaset <- X_ind
  # generalized indexes
  VARS1 <- c("dtaMEI","dtaMHI")
  VARS2 <- c("ddtaMEI","ddtaMHI")
  VARS3 <- c("d2dtaMEI","d2dtaMHI")
  VARS4 <- c(VARS1,VARS2)
  VARS5 <- c(VARS1,VARS3)
  VARS6 <- c(VARS2,VARS3)
  VARS7 <- c(VARS4,VARS3)
  
  VARS <- c("VARS1","VARS2","VARS3","VARS4","VARS5","VARS6","VARS7")
  names <- c("._.MEIMHI",".d.MEIMHI",".d2.MEIMHI","._d.MEIMHI","._d2.MEIMHI",
             ".dd2.MEIMHI","._dd2.MEIMHI")
  
  val_indices <- c()
  clusters <- list()
  times <- data.frame()
  names_method <- c()
  name <- c()
  
  for(i in 1:18){
    #if(abs(det(var(dtaset[,eval(parse(text = VARS[i]))])))>1e-5){
      tryCatch({
        c <- clustInd_svc_aux(dtaset, eval(parse(text = VARS[i])), cluster.method, clus, true_labels)
        name <- paste("svc",names[i],"-",cluster.method,sep="")
        names_method <- c(names_method, name)
        clusters <- c(clusters, list(c$cluster))
        times <- rbind(times, c$time)
        if(!(missing(true_labels))){
          val_indices <- rbind(val_indices,c$valid)
        }
      }, error=function(e){})
    #}
  }
  names(clusters) <- names_method
  row.names(times) <- names_method
  colnames(times)[1] <- "time"
  
  if(missing(true_labels)){
    res <- list("clusters" = clusters, "time" = times) 
  }else{
    row.names(val_indices) <- names_method
    res <- list("val_indices" = as.data.frame(val_indices), "clusters" = clusters, "time" = times) 
  }
  return(res)
}


#########################
### clustInd_spc_aux ####
#########################

clustInd_spc_aux <- function(dtaset, VARS, clus, true_labels){
  
  t0 <- Sys.time()
  met <- specc(na.omit(as.matrix(dtaset[,VARS])), clus)
  clus_part <- met@.Data
  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))
  if(missing(true_labels)){
    res <- list("cluster"=clus_part, "time"=as.numeric(t))
  }else{
    valid <- valid(true_labels,clus_part)
    res <- list("cluster"=clus_part, "valid"=valid, "time"=as.numeric(t))
  }
  return(res)
} 


#####################
### clustInd_spc ####
#####################

## This function applies spectral clustering to a functional data set composed by two bunches 
## of curves and also applies validation techniques to the obtained clustering results.
## Here, execution times are also computed.

clustInd_spc <- function(X_ind,t,rangeval=c(0,1),nbasis=40,norder=4,clus=2,true_labels){
  
  #INPUT
  # X <- functional data object containing the different groups that can be originally observed
  #      (1 to x rows correspond to the observations for the first group, x+1 to y correspond to the second group and so on)
  # t <- sequence for creating the basis
  # rangeval <- a numeric vector of length 2 defining the interval over which the functional data object can be evaluated
  # nbasis <- number of basis functions
  # norder <- order of b-splines
  # clus <- number of clusters
  # true_labels <- data frame containing the correct classification
  
  # dtaset <- ind(X,t,rangeval,nbasis,norder) #multivariate data (applying indexes)
  dtaset <- X_ind
  
  # generalized indexes
  VARS1 <- c("dtaMEI","dtaMHI")
  VARS2 <- c("ddtaMEI","ddtaMHI")
  VARS3 <- c("d2dtaMEI","d2dtaMHI")
  VARS4 <- c(VARS1,VARS2)
  VARS5 <- c(VARS1,VARS3)
  VARS6 <- c(VARS2,VARS3)
  VARS7 <- c(VARS4,VARS3)
  
  VARS <- c("VARS1","VARS2","VARS3","VARS4","VARS5","VARS6","VARS7")
  names <- c("._.MEIMHI",".d.MEIMHI",".d2.MEIMHI","._d.MEIMHI","._d2.MEIMHI",
             ".dd2.MEIMHI","._dd2.MEIMHI")
  
  val_indices <- c()
  clusters <- list()
  times <- data.frame()
  names_method <- c()
  name <- c()
  
  for(i in 1:length(VARS)){
    #if(abs(det(var(dtaset[,eval(parse(text = VARS[i]))])))>1e-5){
      tryCatch({
        c <- clustInd_spc_aux(dtaset, eval(parse(text = VARS[i])), clus, true_labels)
        name <- paste("spc",names[i],sep="")
        names_method <- c(names_method, name)
        clusters <- c(clusters, list(c$cluster))
        times <- rbind(times, c$time)
        if(!(missing(true_labels))){
          val_indices <- rbind(val_indices,c$valid)
        }
      }, error=function(e){})
    #}
  }
  names(clusters) <- names_method
  row.names(times) <- names_method
  colnames(times)[1] <- "time"
  
  if(missing(true_labels)){
    res <- list("clusters" = clusters, "time" = times) 
    
  }else{
    row.names(val_indices) <- names_method
    res <- list("val_indices" = as.data.frame(val_indices), "clusters" = clusters, "time" = times) 
  }
  return(res)
}

##########################################################
################## SIMULATED DATA ########################
##########################################################

## The code for simulating the data is defined in this R script.


sim1 <- function(n=50,t=seq(0,1,length=30),clus2=2){
  X <- matrix(0,n,length(t))
  sigma <- matrix(0,length(t),length(t))
  
  # Building sigma
  for (i in 1:length(t)){
    for (j in 1:length(t)){
      sigma[i,j] <- 0.3*exp((-1/0.3)*abs(t[j]-t[i]));
    }
  }
  
  #centered Gaussian process
  egauss <- mvrnorm(n ,rep(0, length(t)), sigma) #a sample from the specified multivariate normal distribution
  
  sigma2 <- matrix(0,length(t),length(t))
  
  # Construimos sigma
  for (i in 1:length(t)){
    for (j in 1:length(t)){
      sigma2[i,j] <- 0.5*exp((-1/0.2)*abs(t[j]-t[i]));
    }
  }
  
  #centered Gaussian process
  hgauss <- mvrnorm(n ,rep(0, length(t)), sigma) #a sample from the specified multivariate normal distribution
  
  
  ## MODEL 1
  
  X1 <- X
  for (i in 1:n){
    for (j in 1:length(t)){
      X1[i,j] <- 30*t[j]^(3/2)*(1-t[j])+egauss[i,j]
    }
  }
  
  #model clus2
  if(clus2==2){
    X2 <- X
    for (i in 1:n){
      for (j in 1:length(t)){
        X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+0.5+egauss[i,j]
      }
    }
  }
  else if(clus2==3){
    X2 <- X
    for (i in 1:n){
      for (j in 1:length(t)){
        X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+0.75+egauss[i,j]
      }
    }
  }
  else if(clus2==4){
    X2 <- X
    for (i in 1:n){
      for (j in 1:length(t)){
        X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+1.+egauss[i,j]
      }
    }
  }
  else if(clus2==5){
    X2 <- X
    for (i in 1:n){
      for (j in 1:length(t)){
        X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+2*egauss[i,j]
      }
    }
  }
  else if(clus2==6){
    X2 <- X
    for (i in 1:n){
      for (j in 1:length(t)){
        X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+0.25*egauss[i,j]
      }
    }
  }
  else if(clus2==7){
    X2 <- X
    for (i in 1:n){
      for (j in 1:length(t)){
        X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+hgauss[i,j]
      }
    }
  }
  else if(clus2==8){
    X2 <- X
    for (i in 1:n){
      for (j in 1:length(t)){
        X2[i,j] <- 30*t[j]*(1-t[j])+hgauss[i,j]
      }
    }
  }
  else if(clus2==9){
    X2 <- X
    for (i in 1:n){
      for (j in 1:length(t)){
        X2[i,j] <- 30*t[j]*(1-t[j])+egauss[i,j]
      }
    }
  }
  X12 <- rbind(X1, X2)
  return(X12)
}


sim2 <- function(n=50,t=seq(0,1,length.out=150), clus2=11){
  P <- 150 #equidistant points in the interval I=[0,1]
  K <- 100 #components
  m1 <- t*(1-t)
  
  rho <- rep(0, K)
  for(k in 1:K){
    if(k<4)
      rho[k] <- 1/(k+1)
    else
      rho[k] <- 1/(k+1)^2
  }
  
  theta <- matrix(0,K,P)
  for (k in 1:K) {
    if (k%%2 == 0)
      theta[k, ] <- sqrt(2) * sin(k*pi*t)
    else if (k%%2 != 0 && k != 1)
      theta[k, ] <- sqrt( 2 ) * cos(( k-1)*pi* t)
    else
      theta[k, ] <- rep(1, P)
  }
  
  s1 <- 0
  for (k in 1:4) {
    s1 <- s1 + sqrt(rho[k]) * theta[k, ]
  }
  
  s2 <- 0
  for (k in 4:K) {
    s2 <- s2 + sqrt(rho[k]) * theta[k, ]
  }
  
  m2_1 <- m1 + s1
  m2_2 <- m1 + s2
  
  zX <- matrix(0, n, K)
  zY <- matrix(0, n, K)
  
  uX <- matrix(0, n, P)
  uY <- matrix(0, n, P)
  for (i in 1:n) {
    zX[i, ] <- rnorm(K)
    zY[i, ] <- rnorm(K)
    for (k in 1:K) {
      uX[i, ] <- uX[i, ] + sqrt(rho[k]) * (zX[i, k] * theta[k, ])
      uY[i, ] <- uY[i, ] + sqrt(rho[k]) * (zY[i, k] * theta[k, ])
    }
  }
  
  X_ <- matrix(0, n, P) # MODEL 10
  Y1 <- matrix(0, n, P)
  if(clus2==11){
    for (i in 1:n){
      X_[i, ] <- m1 + uX[i, ] 
      Y1[i, ] <- m2_1 + uY[i, ]
    }
  } else if(clus2==12){
    for (i in 1:n){
      X_[i, ] <- m1 + uX[i, ] 
      Y1[i, ] <- m2_2 + uY[i, ]
    }
  }
  X_Y1 <- rbind(X_, Y1)
  return(X_Y1)
}

S1 <-function(n=50,t=seq(0, pi/3, length = 100)){
  #Scenario 1
  a1 <- runif(n,-0.25,0.25)
  b1 <- c(0.3,1,0.2)
  c1 <- c(1/1.3, 1/1.2, 1/4)
  X1 <- matrix(0,n,length(t))
  X2 <- matrix(0,n,length(t))
  X3 <- matrix(0,n,length(t))
  eps <- mvrnorm(n,mu=rep(2,length(t)), Sigma=(0.4^2)*diag(length(t)))
  for (i in 1:n){
    for (j in 1:length(t)){
      X1[i,j] <- a1[i]+b1[1]+c1[1]*sin(1.3*t[j])+(t[j])^3
      X2[i,j] <- a1[i]+b1[2]+c1[2]*sin(1.3*t[j])+(t[j])^3
      X3[i,j] <- a1[i]+b1[3]+c1[3]*sin(1.3*t[j])+(t[j])^3
    }
  }
  
  Y13 <- X1+eps
  Y14 <- X2+eps
  Y15 <- X3+eps
  
  return(rbind(Y13, Y14, Y15))
}


S2 <-function(n=50,t=seq(0, pi/3, length = 100)){
  #Scenario 2
  a1 <- runif(n,-0.5,0.5)
  b1 <- c(1.1, 1.5, 2.2)
  c1 <- c(1.5,1.7,1.9)
  X1 <- matrix(0,n,length(t))
  X2 <- matrix(0,n,length(t))
  X3 <- matrix(0,n,length(t))
  eps <- mvrnorm(n,mu=rep(2,length(t)), Sigma=(0.4^2)*diag(length(t)))
  for (i in 1:n){
    for (j in 1:length(t)){
      X1[i,j] <- a1[i]+b1[1]+sin(c1[1]*pi*t[j])+cos(pi*(t[j])^2)
      X2[i,j] <- a1[i]+b1[2]+sin(c1[2]*pi*t[j])+cos(pi*(t[j])^2)
      X3[i,j] <- a1[i]+b1[3]+sin(c1[3]*pi*t[j])+cos(pi*(t[j])^2)
    }
  }
  
  Y13 <- X1+eps
  Y14 <- X2+eps
  Y15 <- X3+eps
  
  return(rbind(Y13, Y14, Y15))
}


S3 <-function(n=50,t=seq(0, pi/3, length = 100)){
  #Scenario 3
  a1 <- runif(n,-0.25,0.25)
  b1 <- c(1/1.8,1/1.7,1/1.5)
  c1 <- c(1.1,1.4,1.5)
  X1 <- matrix(0,n,length(t))
  X2 <- matrix(0,n,length(t))
  X3 <- matrix(0,n,length(t))
  eps <- mvrnorm(n,mu=rep(2,length(t)), Sigma=(0.3^2)*diag(length(t)))
  for (i in 1:n){
    for (j in 1:length(t)){
      X1[i,j] <- a1[i]+b1[1]*exp(c1[1]*t[j])-(t[j])^3
      X2[i,j] <- a1[i]+b1[2]*exp(c1[2]*t[j])-(t[j])^3
      X3[i,j] <- a1[i]+b1[3]*exp(c1[3]*t[j])-(t[j])^3
    }
  }
  
  Y13 <- X1+eps
  Y14 <- X2+eps
  Y15 <- X3+eps
  
  return(rbind(Y13, Y14, Y15))
}


sim2_mul <- function(n=50,t=seq(0,1,length.out=150), clus2=20){
  P <- 150 #equidistant points in the interval I=[0,1]
  K <- 100 #components
  #X_Y1 <- array(rep(NaN,K*P*2),dim=c(K,P,2)) # datos finales
  X_Y1 <- array(dim=c(K,P,2)) # datos finales
  m1 <- rbind(t*(1-t), 4*t*t*(1-t))
  
  rho <- rep(0, K)
  for(k in 1:K){
    if(k<4)
      rho[k] <- 1/(k+1)
    else
      rho[k] <- 1/(k+1)^2
  }
  
  theta <- matrix(0,K,P)
  for (k in 1:K) {
    if (k%%2 == 0)
      theta[k, ] <- sqrt(2) * sin(k*pi*t)
    else if (k%%2 != 0 && k != 1)
      theta[k, ] <- sqrt( 2 ) * cos(( k-1)*pi* t)
    else
      theta[k, ] <- rep(1, P)
  }
  
  s1 <- 0
  for (k in 1:4) {
    s1 <- s1 + sqrt(rho[k]) * theta[k, ]
  }
  
  s2 <- 0
  for (k in 4:K) {
    s2 <- s2 + sqrt(rho[k]) * theta[k, ]
  }
  
  m2_1 <- rbind(m1[1,] + s1, m1[2,] + s1)
  m2_2 <- rbind(m1[1,] + s2, m1[2,] + s2)
  
  uX1 <- matrix(0, n, P)
  uX2 <- matrix(0, n, P)
  uY1 <- matrix(0, n, P)
  uY2 <- matrix(0, n, P)
  
  mu <- c(0,0) # Mean
  sigma <- matrix(c(1, 0.5, 0.5, 1),
                  2) # Covariance matrix
  
  for (i in 1:n) {
    zX_i <- mvrnorm(K, mu = mu, Sigma = sigma )
    zY_i <- mvrnorm(K, mu = mu, Sigma = sigma )
    for (k in 1:K) {
      uX1[i, ] <- uX1[i, ] + sqrt(rho[k]) * (zX_i[k, 1] * theta[k, ])
      uX2[i, ] <- uX2[i, ] + sqrt(rho[k]) * (zX_i[k, 2] * theta[k, ])
      uY1[i, ] <- uY1[i, ] + sqrt(rho[k]) * (zY_i[k, 1] * theta[k, ])
      uY2[i, ] <- uY2[i, ] + sqrt(rho[k]) * (zY_i[k, 2] * theta[k, ])
    }
  }
  
  X_1 <- matrix(0, n, P) # MODEL 10
  X_2 <- matrix(0, n, P) # MODEL 10
  Y1_1 <- matrix(0, n, P)
  Y1_2 <- matrix(0, n, P)
  if(clus2==20){
    for (i in 1:n){
      X_1[i, ] <- m1[1,] + uX1[i, ] 
      X_2[i, ] <- m1[2,] + uX2[i, ]
      Y1_1[i, ] <- m2_1[1,] + uY1[i, ]
      Y1_2[i, ] <- m2_1[2,] + uY2[i, ]
    }
  } else if(clus2==21){
    for (i in 1:n){
      X_1[i, ] <- m1[1,] + uX1[i, ] 
      X_2[i, ] <- m1[2,] + uX2[i, ]
      Y1_1[i, ] <- m2_2[1,] + uY1[i, ]
      Y1_2[i, ] <- m2_2[2,] + uY2[i, ]
    }
  }
  X_Y1_1 <- rbind(X_1, Y1_1)
  X_Y1_2 <- rbind(X_2, Y1_2)
  
  X_Y1[,,1] <- X_Y1_1
  X_Y1[,,2] <- X_Y1_2
  return(X_Y1)
}


datos_sim <- function(nsim){
  # Data Franco-Pereira and Lillo 2020
  if(nsim>=2 & nsim<=9)
    data <- sim1(clus2 = nsim)
  # Data paper Martino
  else if(nsim==11| nsim==12)
    data <- sim2(clus2 = nsim)
  # Data paper Zambom
  else if(nsim ==13)
    data <- S1()
  else if(nsim ==14)
    data <- S2()
  else if(nsim ==15)
    data <- S3()
  # Multivariate simulations
  else if(nsim ==20|nsim==21)
    data <- sim2_mul(clus2 = nsim)
  return(data)
}


######################################################
########## SIMULATION CODE ###########################
######################################################

create_grid <- function(nsim,t,nbasis,norder,clus,Type){
  #MIRAR SI CONVIENE PONER X EN VEZ DE NSIM PARA PODER USARLA EN CUALQUIER DATASET REAL (pasar X como param)
  X <- datos_sim(nsim)
  clust1 <- clustInd_svc(X,t,c(min(t),max(t)),as.integer(nbasis),cluster.method="kmeans", as.integer(norder),as.integer(clus),true_labels=Type)
  clasif1 <- clust1$val_indices
  t1 <- clust1$time
  clust2 <- clustInd_svc(X,t,c(min(t),max(t)),as.integer(nbasis),cluster.method="mlKmeans", as.integer(norder),as.integer(clus),true_labels=Type)
  clasif2 <- clust2$val_indices
  t2 <- clust2$time
  clust3 <- clustInd_kkmeans(X,t,c(min(t),max(t)),as.integer(nbasis),kernel="rbfdot", as.integer(norder),as.integer(clus),true_labels=Type)
  clasif3 <- clust3$val_indices
  t3 <- clust3$time
  clust4 <- clustInd_kkmeans(X,t,c(min(t),max(t)),as.integer(nbasis),kernel="polydot", as.integer(norder),as.integer(clus),true_labels=Type)
  clasif4 <- clust4$val_indices
  t4 <- clust4$time
  clust5 <- clustInd_spc(X,t,c(min(t),max(t)),as.integer(nbasis), as.integer(norder),as.integer(clus),true_labels=Type)
  clasif5 <- clust5$val_indices
  t5 <- clust5$time
  clust6 <- clustInd_hierarch(X,t,c(min(t),max(t)),as.integer(nbasis),as.integer(norder),"euclidean",as.integer(clus),true_labels=Type)
  clasif6 <- clust6$val_indices
  t6 <- clust6$time
  clust7 <- clustInd_kmeans(X,t,c(min(t),max(t)),as.integer(nbasis),as.integer(norder),"euclidean",as.integer(clus),true_labels=Type)
  clasif7 <- clust7$val_indices
  t7 <- clust7$time
  clust8 <- clustInd_kmeans(X,t,c(min(t),max(t)),as.integer(nbasis),as.integer(norder),"mahalanobis",as.integer(clus),true_labels=Type)
  clasif8 <- clust8$val_indices
  t8 <- clust8$time
  
  clasif1 <- rbind(clasif1,clasif2,clasif3,clasif4,clasif5,clasif6,clasif7,clasif8)
  names <- row.names(clasif1)
  time <- rbind(t1,t2,t3,t4,t5,t6,t7,t8)
  Iter <- rep(1,dim(clasif1)[1])
  clasif <- data.frame(cbind(names,clasif1,time,Iter))
  row.names(clasif)<-c()
  res <- clasif[order(clasif$RI,decreasing = TRUE),]
  return(res)
}

# This function generate the simulation for each generated scenario with all the considered methods

simul_all_fun <- function(nsim,t,nbasis=30,norder=4,clus=2,true_labels,r=100){
  
  #INPUT
  # NSIM <- The number of the consider scenario from function datos_sim
  # t <- sequence for creating the basis
  # nbasis <- number of basis functions
  # norder <- order of b-splines
  # clus <- number of clusters
  # true_labels <- data frame containing the correct classification
  # r <- number of simulations
  
  set.seed(1)
  
  #Here we create reference dataframes to take from them the model names
  
  clasif <- create_grid(nsim,t,nbasis,norder,clus,true_labels)
    
  for(i in 2:r){
    clasif_i <- create_grid(nsim,t,nbasis,norder,clus,true_labels)
    clasif_f <- bind_rows(clasif, clasif_i) %>%
      group_by(names) %>%
      summarise_all(sum)
    clasif <- as.data.frame(clasif_f)
  }  
  #divide Purity, Fmeasure, RI and time by the number of iterations
  clasif[, c(-1,-6)] <- sweep(clasif[, c(-1,-6)], 1, clasif[, 6], "/")
  res <- clasif[order(clasif$RI,decreasing = TRUE),]
  
  return(res)
}

##########################################################
################## EXAMPLES ##############################
##########################################################

# n <- 50
# t <- seq(0,1,length=30)
# Type <- c(rep(1,n), rep(2,n))
# 
# sf3 <- simul_all_fun(3,t,clus=2,true_labels = Type, r=2)
# print(xtable(sf3,digits=5))

# 
# 
# n <- 50
# t <- seq(0,1,length.out=150)
# Type <- c(rep(1,n), rep(2,n))
# 
# sf10 <- simul_all_fun(12,t,clus=2,true_labels = Type, r=100)
# print(xtable(sf10,digits=5))
# 
# 
# n <- 50
# t <- seq(0, pi/3, length = 100)
# Type <- c(rep(1,n), rep(2,n), rep(3,n))
# 
# sf11 <- simul_all_fun(13,t,clus=3,true_labels = Type, r=100)
# print(xtable(sf11,digits=5))


