#==========================================================================================#
# Subroutine for ranking of individual and noise variables main effects 
# by univariate Minimal Depth of a Maximal Subtree (MDMS)
#==========================================================================================#
rsf.main.signif <- function(X,
                            ntree,
                            method,
                            splitrule,
                            importance,
                            B,
                            verbose,
                            seed) {
   
   n <- nrow(X)
   p <- ncol(X) - 2
   X.obs <- X[,-c(1,2)]
   main.obs <- matrix(data=NA, nrow=B, ncol=p)
   main.noise <- matrix(data=NA, nrow=B, ncol=p)
   for (b in 1:B) {
      if (verbose)  cat("Iteration:", b, "\n")
      set.seed(seed=seed[b])
      X.noise <- apply(X=X.obs, MARGIN=2, FUN=function(x) sample(x=x, size=n, replace=FALSE, prob=NULL))
      colnames(X.noise) <- paste(colnames(X.obs), ".noise", sep="")
      Z.obs <- data.frame("time"=X[,1], "event"=X[,2], X.obs)
      Z.noise <- data.frame("time"=X[,1], "event"=X[,2], X.noise)
      set.seed(seed=seed[b])
      bo <- sample(x=1:n, size=n, replace=TRUE, prob=NULL)
      Z.obs.bo <- Z.obs[bo,]
      Z.noise.bo <- Z.noise[bo,]
      
      rsf.obs.bo <- randomForestSRC::rfsrc(formula=survival::Surv(time=time, event=event, type="right") ~ .,
                                           data=Z.obs.bo,
                                           ntree=ntree,
                                           bootstrap="by.root",
                                           mtry=p,
                                           nodesize=3,
                                           splitrule=splitrule,
                                           nsplit=0,
                                           importance=importance,
                                           na.action="na.omit",
                                           proximity=TRUE,
                                           samptype="swr",
                                           forest=TRUE,
                                           var.used="all.trees",
                                           split.depth="all.trees",
                                           membership=TRUE,
                                           statistics=TRUE,
                                           tree.err=TRUE,
                                           seed=seed[b])
      
      rsf.noise.bo <- randomForestSRC::rfsrc(formula=survival::Surv(time=time, event=event, type="right") ~ .,
                                             data=Z.noise.bo,
                                             ntree=ntree,
                                             bootstrap="by.root",
                                             mtry=p,
                                             nodesize=3,
                                             splitrule=splitrule,
                                             nsplit=0,
                                             importance=importance,
                                             na.action="na.omit",
                                             proximity=TRUE,
                                             samptype="swr",
                                             forest=TRUE,
                                             var.used="all.trees",
                                             split.depth="all.trees",
                                             membership=TRUE,
                                             statistics=TRUE,
                                             tree.err=TRUE,
                                             seed=seed[b])
      
      max.obj <- randomForestSRC::max.subtree(object=rsf.obs.bo,
                                              max.order=1,
                                              sub.order=TRUE,
                                              conservative=FALSE)
      main.obs[b,] <- diag(max.obj$sub.order)
      
      max.obj <- randomForestSRC::max.subtree(object=rsf.noise.bo,
                                              max.order=1,
                                              sub.order=TRUE,
                                              conservative=FALSE)
      main.noise[b,] <- diag(max.obj$sub.order)
   }
   colnames(main.obs) <- rsf.obs.bo$xvar.names
   colnames(main.noise) <- rsf.noise.bo$xvar.names
   return(list("boot.obs"=main.obs, "boot.noise"=main.noise))
   
}
#==========================================================================================#



#==========================================================================================#
# Subroutine for ranking of pairwise interactions between individual or noise variables 
# by bivariate interaction Minimal Depth of a Maximal Subtree (IMDMS)
#==========================================================================================#
rsf.int.signif <- function(X,
                           ntree,
                           method,
                           splitrule,
                           importance,
                           B,
                           verbose,
                           seed) {
   
   n <- nrow(X)
   p <- ncol(X) - 2
   X.obs <- X[,-c(1,2)]
   int.mdms.obs <- array(data=NA, dim=c(p, p, B))
   int.mdms.noise <- array(data=NA, dim=c(p, p, B))
   for (b in 1:B) {
      if (verbose)  cat("Iteration:", b, "\n")
      set.seed(seed=seed[b])
      X.noise <- apply(X=X.obs, MARGIN=2, FUN=function(x) sample(x=x, size=n, replace=FALSE, prob=NULL))
      colnames(X.noise) <- paste(colnames(X.obs), ".noise", sep="")
      Z.obs <- data.frame("time"=X[,1], "event"=X[,2], X.obs)
      Z.noise <- data.frame("time"=X[,1], "event"=X[,2], X.noise)
      set.seed(seed=seed[b])
      bo <- sample(x=1:n, size=n, replace=TRUE, prob=NULL)
      Z.obs.bo <- Z.obs[bo,]
      Z.noise.bo <- Z.noise[bo,]
      
      rsf.obs.bo <- randomForestSRC::rfsrc(formula=survival::Surv(time=time, event=event, type="right") ~ .,
                                           data=Z.obs.bo,
                                           ntree=ntree,
                                           bootstrap="by.root",
                                           mtry=p,
                                           nodesize=3,
                                           splitrule=splitrule,
                                           nsplit=0,
                                           importance=importance,
                                           na.action="na.omit",
                                           proximity=TRUE,
                                           samptype="swr",
                                           forest=TRUE,
                                           var.used="all.trees",
                                           split.depth="all.trees",
                                           membership=TRUE,
                                           statistics=TRUE,
                                           tree.err=TRUE,
                                           seed=seed[b])
      
      rsf.noise.bo <- randomForestSRC::rfsrc(formula=survival::Surv(time=time, event=event, type="right") ~ .,
                                             data=Z.noise.bo,
                                             ntree=ntree,
                                             bootstrap="by.root",
                                             mtry=p,
                                             nodesize=3,
                                             splitrule=splitrule,
                                             nsplit=0,
                                             importance=importance,
                                             na.action="na.omit",
                                             proximity=TRUE,
                                             samptype="swr",
                                             forest=TRUE,
                                             var.used="all.trees",
                                             split.depth="all.trees",
                                             membership=TRUE,
                                             statistics=TRUE,
                                             tree.err=TRUE,
                                             seed=seed[b])
      
      int.obs.bo <- randomForestSRC::find.interaction(object=rsf.obs.bo,
                                                      importance=importance,
                                                      method="maxsubtree",
                                                      sorted=TRUE,
                                                      nrep=1,
                                                      na.action="na.omit",
                                                      seed=seed[b],
                                                      verbose=TRUE)
      
      int.noise.bo <- randomForestSRC::find.interaction(object=rsf.noise.bo,
                                                        importance=importance,
                                                        method="maxsubtree",
                                                        sorted=TRUE,
                                                        nrep=1,
                                                        na.action="na.omit",
                                                        seed=seed[b],
                                                        verbose=TRUE)
      
      ord.obs <- order(rownames(int.obs.bo))
      ord.noise <- order(rownames(int.noise.bo))
      int.mdms.obs[,,b] <- int.obs.bo[ord.obs,ord.obs]
      int.mdms.noise[,,b] <- int.noise.bo[ord.noise,ord.noise]
   }
   dimnames(int.mdms.obs) <- list(rownames(int.obs.bo)[ord.obs], colnames(int.obs.bo)[ord.obs], NULL)
   dimnames(int.mdms.noise) <- list(rownames(int.noise.bo)[ord.noise], colnames(int.noise.bo)[ord.noise], NULL)
   return(list("boot.obs"=int.mdms.obs, "boot.noise"=int.mdms.noise))
}
#==========================================================================================#



#==========================================================================================#
#===============#
# Usage         :
#===============#
#                    .onAttach (libname, pkgname)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#==========================================================================================#

.onAttach <- function(libname, pkgname) {
   
   SSver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                     fields="Version")
   packageStartupMessage(paste(pkgname, " ", SSver, sep=""))
   packageStartupMessage("Type IRSF.news() to see new features, changes, and bug fixes")
   
}
#==========================================================================================#
