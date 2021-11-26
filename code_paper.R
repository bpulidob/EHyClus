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
  # nbasis <- number of basis functions --- PARAMETRO QUE NO ESTOY USANDO AHORA
  # norder <- order of b-splines
  
  basisobj <- create.bspline.basis(rangeval=rangeval,norder=norder) #b-spline basis (order 3 polynomials when norder = 4)
  dtaX <- funspline(X,t,basisobj)
  
  dta <- X #original data
  # dta <- dtaX$smooth # smoothed data coefs
  ddta <- dtaX$deriv #first derivatives
  d2dta <- dtaX$deriv2 #second derivatives
  
  # Applying the indexes to the original data
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
  valid <- valid(true_labels, clus_part)
  t <- data.frame(t1-t0)
  res <- list("cluster"=clus_part, "valid"=valid, "time"=as.numeric(t))
  return(res)
}

#########################
### clustInd_hierarch ###
####################[#####

## This function applies hierarchical clustering methods to a data set obtained from applying indexes to
## data composed by a number equal to clus bunches of curves and also applies validation techniques to the 
## obtained clustering results. Execution times are also computed.

clustInd_hierarch <- function(X,t,rangeval=c(0,1),nbasis=40,norder=4,dist="euclidean",clus=2,true_labels){
  
  #INPUT
  # X <- functional data object containing the different groups that can be originally observed
  #      (1 to x rows correspond to the observations for the first group, x+1 to y correspond to the second group and so on)
  # t <- sequence for creating the basis
  # rangeval <- a numeric vector of length 2 defining the interval over which the functional data object can be evaluated
  # nbasis <- number of basis functions
  # norder <- order of b-splines
  # dist <- distance for hierarchical clustering ('euclidean' or 'manhattan')
  # clus <- number of clusters
  # true_labels <- vector containing the correct classification
  
  distances <- c("euclidean","manhattan")
  
  if(is.na(pmatch(dist,distances)))
    stop("invalid distance method")
  
  dtaset <- ind(X,t,rangeval,nbasis,norder) #multivariate dataframe (applying indexes)
  
  # original indexes
  VARS1 <- c("dtaEI","dtaHI")
  VARS2 <- c("ddtaEI","ddtaHI")
  VARS3 <- c("d2dtaEI","d2dtaHI")
  VARS4 <- c(VARS1,VARS2)
  VARS5 <- c(VARS1,VARS3)
  VARS6 <- c(VARS2,VARS3)
  VARS7 <- c(VARS4,VARS3)
  
  # generalized indexes
  VARS8 <- c("dtaMEI","ddtaMEI")
  VARS9 <- c("dtaMEI","d2dtaMEI")
  VARS10 <- c("ddtaMEI","d2dtaMEI")
  VARS11 <- c(VARS8,"d2dtaMEI")
  
  # combining both indexes types
  VARS12 <- c(VARS1,"dtaMEI")
  VARS13 <- c(VARS2,"ddtaMEI")
  VARS14 <- c(VARS3,"d2dtaMEI")
  VARS15 <- c(VARS4,VARS8)
  VARS16 <- c(VARS5,VARS9)
  VARS17 <- c(VARS6,VARS10)
  VARS18 <- c(VARS7,VARS11)
  
  clust_meth <- c("single","complete","average","centroid","ward.D2")
  VARS <- c("VARS1","VARS2","VARS3","VARS4","VARS5","VARS6","VARS7","VARS8","VARS9",
            "VARS10","VARS11","VARS12","VARS13","VARS14","VARS15","VARS16","VARS17","VARS18")
  names <- c("._.EIHI",".d.EIHI",".d2.EIHI","._d.EIHI","._d2.EIHI",".dd2.EIHI","._dd2.EIHI",
             "._d.MEI","._d2.MEI",".dd2.MEI","_dd2.MEI","._.EIHIMEI",".d.EIHIMEI",".d2.EIHIMEI",
             "._d.EIHIMEI","._d2.EIHIMEI",".dd2.EIHIMEI","._dd2.EIHIMEI")
  
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
      if(abs(det(var(dtaset[,eval(parse(text = VARS[i]))])))>1e-5){
        tryCatch({
          c <- clustInd_hierarch_aux(dtaset, eval(parse(text = VARS[i])), clust_meth[m], dist, clus, Type)
          name <- paste(clust_meth[m],names[i],"-",dist,sep="")
          names_method <- c(names_method, name)
          val_indices_i <- rbind(val_indices_i,c$valid)
          clusters_i <- c(clusters_i, list(c$cluster))
          times_i <- rbind(times_i, c$time)
        }, error=function(e){}) 
      }
    }
    row.names(val_indices_i) <- names_method
    names(clusters_i) <- names_method
    row.names(times_i) <- names_method
    colnames(times_i)[1] <- "time"
    
    val_indices <- rbind(val_indices, val_indices_i)
    clusters <- c(clusters, clusters_i)
    times <- rbind(times, times_i)
  }
  res <- list("val_indices" = as.data.frame(val_indices), "clusters" = clusters, "time" = times) 
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
  valid <- valid(true_labels,clus_part)
  t <- data.frame(t1-t0)
  res <- list("cluster"=clus_part, "valid"=valid, "time"=as.numeric(t))
  return(res)
}  

#######################
### clustInd_kmeans ###
#######################

## This function applies partitional clustering technique (k-means) to a data set obtained from applying indexes to
## data composed by two bunches of curves and also applies validation techniques to the obtained clustering results.
## Here, execution times are also computed.

clustInd_kmeans <- function(X,t,rangeval=c(0,1),nbasis=40,norder=4,dist='euclidean',clus=2,true_labels){
  
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
  
  dtaset <- ind(X,t,rangeval,nbasis,norder) #multivariate data (applying indexes)
  
  # original indexes
  VARS1 <- c("dtaEI","dtaHI")
  VARS2 <- c("ddtaEI","ddtaHI")
  VARS3 <- c("d2dtaEI","d2dtaHI")
  VARS4 <- c(VARS1,VARS2)
  VARS5 <- c(VARS1,VARS3)
  VARS6 <- c(VARS2,VARS3)
  VARS7 <- c(VARS4,VARS3)
  
  # generalized indexes
  VARS8 <- c("dtaMEI","ddtaMEI")
  VARS9 <- c("dtaMEI","d2dtaMEI")
  VARS10 <- c("ddtaMEI","d2dtaMEI")
  VARS11 <- c(VARS8,"d2dtaMEI")
  
  # combining both indexes types
  VARS12 <- c(VARS1,"dtaMEI")
  VARS13 <- c(VARS2,"ddtaMEI")
  VARS14 <- c(VARS3,"d2dtaMEI")
  VARS15 <- c(VARS4,VARS8)
  VARS16 <- c(VARS5,VARS9)
  VARS17 <- c(VARS6,VARS10)
  VARS18 <- c(VARS7,VARS11)
  
  VARS <- c("VARS1","VARS2","VARS3","VARS4","VARS5","VARS6","VARS7","VARS8","VARS9",
            "VARS10","VARS11","VARS12","VARS13","VARS14","VARS15","VARS16","VARS17","VARS18")
  names <- c("._.EIHI",".d.EIHI",".d2.EIHI","._d.EIHI","._d2.EIHI",".dd2.EIHI","._dd2.EIHI",
             "._d.MEI","._d2.MEI",".dd2.MEI","_dd2.MEI","._.EIHIMEI",".d.EIHIMEI",".d2.EIHIMEI",
             "._d.EIHIMEI","._d2.EIHIMEI",".dd2.EIHIMEI","._dd2.EIHIMEI")
  
  val_indices <- c()
  clusters <- list()
  times <- data.frame()
  names_method <- c()
  name <- c()
  
  for(i in 1:18){
    if(abs(det(var(dtaset[,eval(parse(text = VARS[i]))])))>1e-5){
      tryCatch({
        c <- clustInd_kmeans_aux(dtaset, eval(parse(text = VARS[i])), dist, clus, Type)
        name <- paste("kmeans",names[i],"-",dist,sep="")
        names_method <- c(names_method, name)
        val_indices <- rbind(val_indices,c$valid)
        clusters <- c(clusters, list(c$cluster))
        times <- rbind(times, c$time)
      }, error=function(e){}) 
    }
  }
  row.names(val_indices) <- names_method
  names(clusters) <- names_method
  row.names(times) <- names_method
  colnames(times)[1] <- "time"
  
  res <- list("val_indices" = as.data.frame(val_indices), "clusters" = clusters, "time" = times) 
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
  valid <- valid(true_labels,clus_part)
  t <- data.frame(t1-t0)
  res <- list("cluster"=clus_part, "valid"=valid, "time"=as.numeric(t))
  return(res)
} 

##########################
### clustInd_kkmeans ####
#########################

## This function applies kernel k-means to a functional data set composed by two bunches 
## of curves and also applies validation techniques to the obtained clustering results.
## Here, execution times are also computed.

clustInd_kkmeans <- function(X,t,rangeval=c(0,1),nbasis=40,kernel ="rbfdot",norder=4,clus=2,true_labels){
  
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
  
  dtaset <- ind(X,t,rangeval,nbasis,norder) #multivariate data (applying indexes)
  
  # original indexes
  VARS1 <- c("dtaEI","dtaHI")
  VARS2 <- c("ddtaEI","ddtaHI")
  VARS3 <- c("d2dtaEI","d2dtaHI")
  VARS4 <- c(VARS1,VARS2)
  VARS5 <- c(VARS1,VARS3)
  VARS6 <- c(VARS2,VARS3)
  VARS7 <- c(VARS4,VARS3)
  
  # generalized indexes
  VARS8 <- c("dtaMEI","ddtaMEI")
  VARS9 <- c("dtaMEI","d2dtaMEI")
  VARS10 <- c("ddtaMEI","d2dtaMEI")
  VARS11 <- c(VARS8,"d2dtaMEI")
  
  # combining both indexes types
  VARS12 <- c(VARS1,"dtaMEI")
  VARS13 <- c(VARS2,"ddtaMEI")
  VARS14 <- c(VARS3,"d2dtaMEI")
  VARS15 <- c(VARS4,VARS8)
  VARS16 <- c(VARS5,VARS9)
  VARS17 <- c(VARS6,VARS10)
  VARS18 <- c(VARS7,VARS11)
  
  VARS <- c("VARS1","VARS2","VARS3","VARS4","VARS5","VARS6","VARS7","VARS8","VARS9",
            "VARS10","VARS11","VARS12","VARS13","VARS14","VARS15","VARS16","VARS17","VARS18")
  names <- c("._.EIHI",".d.EIHI",".d2.EIHI","._d.EIHI","._d2.EIHI",".dd2.EIHI","._dd2.EIHI",
             "._d.MEI","._d2.MEI",".dd2.MEI","_dd2.MEI","._.EIHIMEI",".d.EIHIMEI",".d2.EIHIMEI",
             "._d.EIHIMEI","._d2.EIHIMEI",".dd2.EIHIMEI","._dd2.EIHIMEI")
  
  val_indices <- c()
  clusters <- list()
  times <- data.frame()
  names_method <- c()
  name <- c()
  
  for(i in 1:18){
    if(abs(det(var(dtaset[,eval(parse(text = VARS[i]))])))>1e-5){
      tryCatch({
        c <- clustInd_kkmeans_aux(dtaset, eval(parse(text = VARS[i])), kernel, clus, Type)
        name <- paste("kkmeans",names[i],"-",kernel,sep="")
        names_method <- c(names_method, name)
        val_indices <- rbind(val_indices,c$valid)
        clusters <- c(clusters, list(c$cluster))
        times <- rbind(times, c$time)
      }, error=function(e){})  
    }
  }
  row.names(val_indices) <- names_method
  names(clusters) <- names_method
  row.names(times) <- names_method
  colnames(times)[1] <- "time"
  
  res <- list("val_indices" = as.data.frame(val_indices), "clusters" = clusters, "time" = times) 
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
  valid <- valid(true_labels,clus_part)
  t <- data.frame(t1-t0)
  res <- list("cluster"=clus_part, "valid"=valid, "time"=as.numeric(t))
  return(res)
} 

#####################
### clustInd_svc ####
#####################

## This function applies support vector clustering to a functional data set composed by two bunches 
## of curves and also applies validation techniques to the obtained clustering results.
## Here, execution times are also computed.

clustInd_svc <- function(X,t,rangeval=c(0,1),nbasis=40,cluster.method="kmeans",norder=4,clus=2,true_labels){
  
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
  
  dtaset <- ind(X,t,rangeval,nbasis,norder) #multivariate data (applying indexes)
  
  # original indexes
  VARS1 <- c("dtaEI","dtaHI")
  VARS2 <- c("ddtaEI","ddtaHI")
  VARS3 <- c("d2dtaEI","d2dtaHI")
  VARS4 <- c(VARS1,VARS2)
  VARS5 <- c(VARS1,VARS3)
  VARS6 <- c(VARS2,VARS3)
  VARS7 <- c(VARS4,VARS3)
  
  # generalized indexes
  VARS8 <- c("dtaMEI","ddtaMEI")
  VARS9 <- c("dtaMEI","d2dtaMEI")
  VARS10 <- c("ddtaMEI","d2dtaMEI")
  VARS11 <- c(VARS8,"d2dtaMEI")
  
  # combining both indexes types
  VARS12 <- c(VARS1,"dtaMEI")
  VARS13 <- c(VARS2,"ddtaMEI")
  VARS14 <- c(VARS3,"d2dtaMEI")
  VARS15 <- c(VARS4,VARS8)
  VARS16 <- c(VARS5,VARS9)
  VARS17 <- c(VARS6,VARS10)
  VARS18 <- c(VARS7,VARS11)
  
  VARS <- c("VARS1","VARS2","VARS3","VARS4","VARS5","VARS6","VARS7","VARS8","VARS9",
            "VARS10","VARS11","VARS12","VARS13","VARS14","VARS15","VARS16","VARS17","VARS18")
  names <- c("._.EIHI",".d.EIHI",".d2.EIHI","._d.EIHI","._d2.EIHI",".dd2.EIHI","._dd2.EIHI",
             "._d.MEI","._d2.MEI",".dd2.MEI","_dd2.MEI","._.EIHIMEI",".d.EIHIMEI",".d2.EIHIMEI",
             "._d.EIHIMEI","._d2.EIHIMEI",".dd2.EIHIMEI","._dd2.EIHIMEI")
  
  val_indices <- c()
  clusters <- list()
  times <- data.frame()
  names_method <- c()
  name <- c()
  
  for(i in 1:18){
    if(abs(det(var(dtaset[,eval(parse(text = VARS[i]))])))>1e-5){
      tryCatch({
        c <- clustInd_svc_aux(dtaset, eval(parse(text = VARS[i])), cluster.method, clus, Type)
        name <- paste("svc",names[i],"-",cluster.method,sep="")
        names_method <- c(names_method, name)
        val_indices <- rbind(val_indices,c$valid)
        clusters <- c(clusters, list(c$cluster))
        times <- rbind(times, c$time)
      }, error=function(e){})
    }
  }
  row.names(val_indices) <- names_method
  names(clusters) <- names_method
  row.names(times) <- names_method
  colnames(times)[1] <- "time"
  
  res <- list("val_indices" = as.data.frame(val_indices), "clusters" = clusters, "time" = times) 
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
  valid <- valid(true_labels,clus_part)
  t <- data.frame(t1-t0)
  res <- list("cluster"=clus_part, "valid"=valid, "time"=as.numeric(t))
  return(res)
} 


#####################
### clustInd_spc ####
#####################

## This function applies spectral clustering to a functional data set composed by two bunches 
## of curves and also applies validation techniques to the obtained clustering results.
## Here, execution times are also computed.

clustInd_spc <- function(X,t,rangeval=c(0,1),nbasis=40,norder=4,clus=2,true_labels){
  
  #INPUT
  # X <- functional data object containing the different groups that can be originally observed
  #      (1 to x rows correspond to the observations for the first group, x+1 to y correspond to the second group and so on)
  # t <- sequence for creating the basis
  # rangeval <- a numeric vector of length 2 defining the interval over which the functional data object can be evaluated
  # nbasis <- number of basis functions
  # norder <- order of b-splines
  # clus <- number of clusters
  # true_labels <- data frame containing the correct classification
  
  dtaset <- ind(X,t,rangeval,nbasis,norder) #multivariate data (applying indexes)
  
  # original indexes
  VARS1 <- c("dtaEI","dtaHI")
  VARS2 <- c("ddtaEI","ddtaHI")
  VARS3 <- c("d2dtaEI","d2dtaHI")
  VARS4 <- c(VARS1,VARS2)
  VARS5 <- c(VARS1,VARS3)
  VARS6 <- c(VARS2,VARS3)
  VARS7 <- c(VARS4,VARS3)
  
  # generalized indexes
  VARS8 <- c("dtaMEI","ddtaMEI")
  VARS9 <- c("dtaMEI","d2dtaMEI")
  VARS10 <- c("ddtaMEI","d2dtaMEI")
  VARS11 <- c(VARS8,"d2dtaMEI")
  
  # combining both indexes types
  VARS12 <- c(VARS1,"dtaMEI")
  VARS13 <- c(VARS2,"ddtaMEI")
  VARS14 <- c(VARS3,"d2dtaMEI")
  VARS15 <- c(VARS4,VARS8)
  VARS16 <- c(VARS5,VARS9)
  VARS17 <- c(VARS6,VARS10)
  VARS18 <- c(VARS7,VARS11)
  
  VARS <- c("VARS1","VARS2","VARS3","VARS4","VARS5","VARS6","VARS7","VARS8","VARS9",
            "VARS10","VARS11","VARS12","VARS13","VARS14","VARS15","VARS16","VARS17","VARS18")
  names <- c("._.EIHI",".d.EIHI",".d2.EIHI","._d.EIHI","._d2.EIHI",".dd2.EIHI","._dd2.EIHI",
             "._d.MEI","._d2.MEI",".dd2.MEI","_dd2.MEI","._.EIHIMEI",".d.EIHIMEI",".d2.EIHIMEI",
             "._d.EIHIMEI","._d2.EIHIMEI",".dd2.EIHIMEI","._dd2.EIHIMEI")
  
  val_indices <- c()
  clusters <- list()
  times <- data.frame()
  names_method <- c()
  name <- c()
  
  for(i in 1:length(VARS)){
    if(abs(det(var(dtaset[,eval(parse(text = VARS[i]))])))>1e-5){
      tryCatch({
        c <- clustInd_spc_aux(dtaset, eval(parse(text = VARS[i])), clus, Type)
        name <- paste("spc",names[i],sep="")
        names_method <- c(names_method, name)
        val_indices <- rbind(val_indices,c$valid)
        clusters <- c(clusters, list(c$cluster))
        times <- rbind(times, c$time)
      }, error=function(e){})
    }
  }
  row.names(val_indices) <- names_method
  names(clusters) <- names_method
  row.names(times) <- names_method
  colnames(times)[1] <- "time"
  
  res <- list("val_indices" = as.data.frame(val_indices), "clusters" = clusters, "time" = times) 
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


datos_sim <- function(nsim){
  if(nsim>=2 & nsim<=9)
    data <- sim1(clus2 = nsim)
  else if(nsim==11| nsim==12)
    data <- sim2(clus2 = nsim)
  else if(nsim ==13)
    data <- S1()
  else if(nsim ==14)
    data <- S2()
  else if(nsim ==15)
    data <- S3()
  return(data)
}

######################################################
########## SIMULATION CODE ###########################
######################################################

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
  

  #Here we create reference dataframes to take from them the model names
  X <- datos_sim(nsim)
  ind_ref1 <- clustInd_svc(X,t,c(min(t),max(t)),nbasis,cluster.method="kmeans",norder,clus,true_labels)$val_indices
  d1 <- dim(ind_ref1)
  ind_ref2 <- clustInd_svc(X,t,c(min(t),max(t)),nbasis,cluster.method="kernKmeans",norder,clus,true_labels)$val_indices
  d2 <- dim(ind_ref2)
  ind_ref3 <- clustInd_kkmeans(X,t,c(min(t),max(t)),nbasis,kernel="rbfdot",norder,clus,true_labels)$val_indices
  d3 <- dim(ind_ref3)
  ind_ref4 <- clustInd_kkmeans(X,t,c(min(t),max(t)),nbasis,kernel="polydot",norder,clus,true_labels)$val_indices
  d4 <- dim(ind_ref4)
  ind_ref5 <- clustInd_spc(X,t,c(min(t),max(t)),nbasis,norder,clus,true_labels)$val_indices
  d5 <- dim(ind_ref5)
  ind_ref6 <- clustInd_hierarch(X,t,c(min(t),max(t)),nbasis,norder,"euclidean",clus,true_labels)$val_indices
  d6 <- dim(ind_ref6)
  ind_ref7 <- clustInd_kmeans(X,t,c(min(t),max(t)),nbasis,norder,"euclidean",clus,true_labels)$val_indices
  d7 <- dim(ind_ref7)
  ind_ref8 <- clustInd_kmeans(X,t,c(min(t),max(t)),nbasis,norder,"mahalanobis",clus,true_labels)$val_indices
  d8 <- dim(ind_ref8)
  ind_ref <- rbind(ind_ref1,ind_ref2,ind_ref3,ind_ref4,ind_ref5,ind_ref6,ind_ref7,ind_ref8)
  nr <- dim(ind_ref)[1] # number of rows
  rn <- row.names(ind_ref) # row names
  
  ### store mean results
  
  Purity <- rep(0,nr)
  Fmeasure <- rep(0,nr)
  RI <- rep(0,nr)
  time_t <- rep(0,nr)
  
  df1 <- data.frame(Purity,Fmeasure,RI,row.names = rn)
  df2 <- data.frame(time_t,row.names = rn)
  
  cont_t <- 0 #number of simulations carried out (initially set to 0)
  cont <- 0 #number of simulations finally considered (initially set to 0)
  
  for (s in 1:(2*r)){
    while(cont<r){
      X <- datos_sim(nsim)
      cont_t <- cont_t+1
      
      clust1 <- clustInd_svc(X,t,c(min(t),max(t)),as.integer(nbasis),cluster.method="kmeans", as.integer(norder),as.integer(clus),true_labels=Type)
      clasif1 <- clust1$val_indices
      t1 <- clust1$time
      clust2 <- clustInd_svc(X,t,c(min(t),max(t)),as.integer(nbasis),cluster.method="kernkmeans", as.integer(norder),as.integer(clus),true_labels=Type)
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
      
      if((d1==dim(clasif1)) && (d2==dim(clasif2)) && (d3==dim(clasif3)) && (d4==dim(clasif4))  
         && (d5==dim(clasif5)) && (d6==dim(clasif6)) && (d7==dim(clasif7))  && (d8==dim(clasif8))){ #discarding simulated data obtaining different set of variables than the first one
        clasif <- rbind(clasif1,clasif2,clasif3,clasif4,clasif5,clasif6,clasif7,clasif8)
        tt <- rbind(t1,t2,t3,t4,t5,t6,t7,t8)
        df1 <- df1+clasif
        df2 <- df2+tt
        cont <- cont+1
      }
    }
  }  
  res <- list("models_names"= rn,"mean_coef"=(df1[order(df1$RI,decreasing = TRUE),])/cont,
              "mean_time"=df2/cont_t, "num_sim"=cont)
  
  return(res)
}

##########################################################
################## EXAMPLES ##############################
##########################################################

# n <- 50
# t <- seq(0,1,length=30)
# Type <- c(rep(1,n), rep(2,n))
# 
# sf3 <- simul_all_fun(4,t,clus=2,true_labels = Type, r=100)
# print(xtable(sf3$mean_coef,digits=3))
# print(xtable(sf3$mean_time,digits=4))
# m3 <- merge(sf3$mean_coef, sf3$mean_time, by=0)
# m3_ord <- m3[order(m3$RI,decreasing = TRUE),]
# print(xtable(m3_ord,digits=5))
# 
# 
# n <- 50
# t <- seq(0,1,length.out=150)
# Type <- c(rep(1,n), rep(2,n))
# 
# sf10 <- simul_all_fun(12,t,clus=2,true_labels = Type, r=100)
# print(xtable(sf10$mean_coef,digits=4))
# print(xtable(sf10$mean_time,digits=4))
# m10 <- merge(sf10$mean_coef, sf10$mean_time, by=0)
# m10_ord <- m10[order(m10$RI,decreasing = TRUE),]
# print(xtable(m10_ord,digits=5))
# 
# 
# n <- 50
# t <- seq(0, pi/3, length = 100)
# Type <- c(rep(1,n), rep(2,n), rep(3,n))
# 
# sf11 <- simul_all_fun(13,t,clus=3,true_labels = Type, r=100)
# print(xtable(sf11$mean_coef,digits=4))
# print(xtable(sf11mean_time,digits=4))
# m11 <- merge(sf11$mean_coef, sf11$mean_time, by=0)
# m11_ord <- m11[order(m11$RI,decreasing = TRUE),]
# print(xtable(m11_ord,digits=5))
