#==========================================================================================#
# Ranking of individual and noise variables main effects 
# by univariate Minimal Depth of a Maximal Subtree (MDMS)
#==========================================================================================#
rsf.main <- function(X,
                     ntree = 1000,
                     method = "mdms",
                     splitrule = "logrank",
                     importance = "random",
                     B,
                     ci = 90,
                     parallel = FALSE,
                     conf = NULL,
                     verbose = TRUE,
                     seed = NULL) {
   
   if (!parallel) {
      if (is.null(seed)) {
         digits <- getOption("digits")
         seed <- runif(n=B, min=1, max=2) * 10^(digits-2)
      } else {
         seed <- (0:(B-1)) + seed
      }
      rsf.obj <- rsf.main.signif(X=X,
                                 ntree=ntree,
                                 method=method,
                                 splitrule=splitrule,
                                 importance=importance,
                                 B=B,
                                 verbose=verbose,
                                 seed=seed)
      rsf.obs.boot <- rsf.obj$boot.obs
      rsf.noise.boot <- rsf.obj$boot.noise
   } else {
      if (conf$type == "SOCKET") {
         cl <- parallel::makeCluster(spec=conf$spec,
                                     type="PSOCK",
                                     homogeneous=conf$homo,
                                     outfile=conf$outfile,
                                     verbose=conf$verbose)
         cpus <- length(conf$spec)
      } else if (conf$type == "MPI") {
         cl <- parallel::makeCluster(spec=conf$spec,
                                     type="MPI",
                                     homogeneous=conf$homo,
                                     outfile=conf$outfile,
                                     verbose=conf$verbose)
         cpus <- conf$spec
      } else {
         stop("Unrecognized cluster type\n")
      }
      parallel::clusterSetRNGStream(cl=cl, iseed=seed)
      parallel::clusterEvalQ(cl=cl, expr=library("randomForestSRC"))
      parallel::clusterEvalQ(cl=cl, expr=library("survival"))
      parallel::clusterExport(cl=cl,
                              varlist=c("rsf.main.signif"),
                              envir=.GlobalEnv)
      rsf.cl <- parallel::clusterCall(cl=cl, fun=rsf.main.signif,
                                      X=X,
                                      ntree=ntree,
                                      method=method,
                                      splitrule=splitrule,
                                      importance=importance,
                                      B=ceiling(B/cpus),
                                      verbose=verbose,
                                      seed=NULL)
      parallel::stopCluster(cl)
      rsf.obs.boot <- matrix(data=NA, nrow=0, ncol=ncol(rsf.cl[[1]]$boot.obs))
      rsf.noise.boot <- matrix(data=NA, nrow=0, ncol=ncol(rsf.cl[[1]]$boot.noise))
      for (b in 1:cpus) {
         rsf.obs.boot <-  rbind(rsf.obs.boot, rsf.cl[[b]]$boot.obs)
         rsf.noise.boot <-  rbind(rsf.noise.boot, rsf.cl[[b]]$boot.noise)
      }
   }
   theta <- (100-ci)/200
   ranks.obs.mean <- apply(X=rsf.obs.boot, MARGIN=2, FUN=mean)
   ranks.noise.mean <- apply(X=rsf.noise.boot, MARGIN=2, FUN=mean)
   ranks.obs.se <- apply(X=rsf.obs.boot, MARGIN=2, FUN=sd)
   ranks.noise.se <- apply(X=rsf.noise.boot, MARGIN=2, FUN=sd)
   ranks.obs.bpci <- apply(X=rsf.obs.boot, MARGIN=2, FUN=function(x) quantile(x=x, probs=c(theta, 1-theta)))
   ranks.noise.bpci <- apply(X=rsf.noise.boot, MARGIN=2, FUN=function(x) quantile(x=x, probs=c(theta, 1-theta)))
   names(ranks.noise.mean) <- gsub(pattern=".noise", replacement="", x=names(ranks.noise.mean), ignore.case=F, fixed=F)
   w <- pmatch(x=names(ranks.obs.mean), table=names(ranks.noise.mean))
   ranks.noise.mean <- ranks.noise.mean[w]
   ranks.noise.se <- ranks.noise.se[w]
   ranks.noise.bpci <- ranks.noise.bpci[,w]
   if (method == "mdms") {
      mat.ranks <- data.frame("obs.mean"=ranks.obs.mean, "obs.se"=ranks.obs.se,
                              "obs.LBCI"=ranks.obs.bpci[1,], "obs.UBCI"=ranks.obs.bpci[2,],
                              "noise.mean"=ranks.noise.mean, "noise.se"=ranks.noise.se,
                              "noise.LBCI"=ranks.noise.bpci[1,], "noise.UBCI"=ranks.noise.bpci[2,],
                              "signif.1SE"=(ranks.obs.mean + ranks.obs.se < ranks.noise.mean),
                              "signif.CI"=(ranks.obs.bpci[2,] < ranks.noise.bpci[1,]))
      ord <- order(mat.ranks[,"obs.mean"], decreasing=F)
   } else {
      stop("Unmatched method \n")
   }
   mat.ranks <- mat.ranks[ord,]
   return(mat.ranks)
   
}
#==========================================================================================#



#==========================================================================================#
# Ranking of pairwise interactions between individual or noise variables 
# by bivariate interaction Minimal Depth of a Maximal Subtree (IMDMS)
#==========================================================================================#
rsf.int <- function(X,
                    ntree = 1000,
                    method = "imdms",
                    splitrule = "logrank",
                    importance = "random",
                    B,
                    ci = 90,
                    parallel = FALSE,
                    conf = NULL,
                    verbose = TRUE,
                    seed = NULL) {
   
   p <- ncol(X) - 2
   if (!parallel) {
      if (is.null(seed)) {
         digits <- getOption("digits")
         seed <- runif(n=B, min=1, max=2) * 10^(digits-2)
      } else {
         seed <- (0:(B-1)) + seed
      }
      rsf.obj <- rsf.int.signif(X=X,
                                ntree=ntree,
                                method=method,
                                splitrule=splitrule,
                                importance=importance,
                                B=B,
                                verbose=verbose,
                                seed=seed)
      rsf.obs.boot <- rsf.obj$boot.obs
      rsf.noise.boot <- rsf.obj$boot.noise
   } else {
      if (conf$type == "SOCKET") {
         cl <- parallel::makeCluster(spec=conf$spec,
                                     type="PSOCK",
                                     homogeneous=conf$homo,
                                     outfile=conf$outfile,
                                     verbose=conf$verbose)
         cpus <- length(conf$spec)
      } else if (conf$type == "MPI") {
         cl <- parallel::makeCluster(spec=conf$spec,
                                     type="MPI",
                                     homogeneous=conf$homo,
                                     outfile=conf$outfile,
                                     verbose=conf$verbose)
         cpus <- conf$spec
      } else {
         stop("Unrecognized cluster type\n")
      }
      parallel::clusterSetRNGStream(cl=cl, iseed=seed)
      parallel::clusterEvalQ(cl=cl, expr=library("randomForestSRC"))
      parallel::clusterEvalQ(cl=cl, expr=library("survival"))
      parallel::clusterExport(cl=cl,
                              varlist=c("rsf.int.signif"),
                              envir=.GlobalEnv)
      rsf.cl <- parallel::clusterCall(cl=cl, fun=rsf.int.signif,
                                      X=X,
                                      ntree=ntree,
                                      method=method,
                                      splitrule=splitrule,
                                      importance=importance,
                                      B=ceiling(B/cpus),
                                      verbose=verbose,
                                      seed=NULL)
      parallel::stopCluster(cl)
      if (method == "imdms") {
         rsf.obs.boot <- array(data=NA, dim=c(dim(rsf.cl[[1]]$boot.obs)[1], dim(rsf.cl[[1]]$boot.obs)[2], 0))
         rsf.noise.boot <- array(data=NA, dim=c(dim(rsf.cl[[1]]$boot.noise)[1], dim(rsf.cl[[1]]$boot.noise)[2], 0))
         for (b in 1:cpus) {
            rsf.obs.boot <-  abind::abind(rsf.obs.boot, rsf.cl[[b]]$boot.obs)
            rsf.noise.boot <-  abind::abind(rsf.noise.boot, rsf.cl[[b]]$boot.noise)
         }
      } else {
         stop("Unmatched method \n")
      }
   }
   theta <- (100-ci)/200
   if (method == "imdms") {
      int.obs.mean <- apply(X=rsf.obs.boot, MARGIN=1:2, FUN=mean)
      int.noise.mean <- apply(X=rsf.noise.boot, MARGIN=1:2, FUN=mean)
      int.obs.se <- apply(X=rsf.obs.boot, MARGIN=1:2, FUN=sd)
      int.noise.se <- apply(X=rsf.noise.boot, MARGIN=1:2, FUN=sd)
      int.obs.bpci <- apply(X=rsf.obs.boot, MARGIN=1:2, FUN=function(x) quantile(x=x, probs=c(theta, 1-theta)))
      int.noise.bpci <- apply(X=rsf.noise.boot, MARGIN=1:2, FUN=function(x) quantile(x=x, probs=c(theta, 1-theta)))
      int.obs.bpci <- aperm(a=int.obs.bpci, perm=c(2,3,1))
      int.noise.bpci <- aperm(a=int.noise.bpci, perm=c(2,3,1))
      vo.mean <- numeric(choose(n=p, k=2))
      vn.mean <- numeric(choose(n=p, k=2))
      vo.se <- numeric(choose(n=p, k=2))
      vn.se <- numeric(choose(n=p, k=2))
      mo.bpci <- matrix(data=NA, nrow=4, ncol=choose(n=p, k=2))
      mn.bpci <- matrix(data=NA, nrow=4, ncol=choose(n=p, k=2))
      k <- 1
      for (i in 2:p) {
         for (j in 1:(i-1)) {
            vmean <- c(int.obs.mean[i,j], int.obs.mean[j,i])
            vse <- c(int.obs.se[i,j], int.obs.se[j,i])
            mbpci <- cbind(int.obs.bpci[i,j,], int.obs.bpci[j,i,])
            vo.mean[k] <- vmean[which.min(vmean)]
            vo.se[k] <- vse[which.min(vmean)]
            mo.bpci[,k] <- mbpci[,which.min(vmean)]
            names(vo.mean)[k] <- paste(colnames(int.obs.mean)[j], rownames(int.obs.mean)[i], sep=":")
            vmean <- c(int.noise.mean[i,j], int.noise.mean[j,i])
            vse <- c(int.noise.se[i,j], int.noise.se[j,i])
            mbpci <- cbind(int.noise.bpci[i,j,], int.noise.bpci[j,i,])
            vn.mean[k] <- vmean[which.min(vmean)]
            vn.se[k] <- vse[which.min(vmean)]
            mn.bpci[,k] <- mbpci[,which.min(vmean)]
            names(vn.mean)[k] <- paste(colnames(int.noise.mean)[j], rownames(int.noise.mean)[i], sep=":")
            k <- k+1
         }
      }
      mat.int <- data.frame("obs.mean"=vo.mean, "obs.se"=vo.se,
                            "obs.LBCI"=mo.bpci[1,], "obs.UBCI"=mo.bpci[2,],
                            "noise.mean"=vn.mean, "noise.se"=vn.se,
                            "noise.LBCI"=mn.bpci[1,], "noise.UBCI"=mn.bpci[2,],
                            "signif.1SE"=(vo.mean + vo.se < vn.mean),
                            "signif.CI"=(mo.bpci[2,] < mn.bpci[1,]))
      ord <- order(mat.int[,"obs.mean"], decreasing=F)
   } else {
      stop("Unmatched method \n")
   }
   mat.int <- mat.int[ord,]
   return(mat.int)
   
}
#==========================================================================================#



#==========================================================================================#
# Fits a Proportional Hazards Time-To-Event Regression Model saturated with first order terms
# Computes p-values of significance of regression coefficients of main effects in a Cox-PH model
#==========================================================================================#
cph.main <- function (X, main.term) {
   
   p <- ncol(X) - 2
   fmla.main <- as.formula(paste("survival::Surv(time=", colnames(X)[1], ", event=", colnames(X)[2], ", type=\"right\") ~ .", sep=""))
   
   P.cph.main <- numeric(p)
   names(P.cph.main) <- main.term
   for (j in 1:p) {
      Z <- X[,c(colnames(X)[1], colnames(X)[2], main.term[j])]
      coxfit <- tryCatch({survival::coxph(fmla.main, data=Z, model=T, x=T, y=T)}, error=function(w){NULL}, warning=function(w){NULL})
      if (is.null(coxfit)) {
         P.cph.main[j] <- NA
      } else {
         P.cph.main[j] <- summary(coxfit)$coefficients[1,"Pr(>|z|)",drop=T]
      }
   }
   P.cph.bh.main <- P.cph.main
   nna <- !is.na(P.cph.main)
   P.cph.bh.main[nna] <- p.adjust(p=P.cph.main[nna], method="BH")
   P.cph.bh.main[!nna] <- NA
   
   return(list("raw"=P.cph.main, "fdr"=P.cph.bh.main))
}
#==========================================================================================#



#==========================================================================================#
# Fits a Proportional Hazards Time-To-Event Regression Model saturated with first and second order terms
# Computes p-values of significance of regression coefficients of pairwise interaction effects in a Cox-PH model
#==========================================================================================#
cph.int <- function (X, int.term) {
   
   p <- ncol(X) - 2
   int.term.list <- strsplit(x=int.term, split=":")
   fmla.int <- as.formula(paste("survival::Surv(time=", colnames(X)[1], ", event=", colnames(X)[2], ", type=\"right\") ~ . + (.)^2", sep=""))
   
   P.cph.int <- numeric(p*(p-1)/2)
   names(P.cph.int) <- int.term
   for (j in 1:(p*(p-1)/2)) {
      Z <- X[,c(colnames(X)[1], colnames(X)[2], int.term.list[[j]])]
      coxfit <- tryCatch({survival::coxph(fmla.int, data=Z, model=T, x=T, y=T)}, error=function(w){NULL}, warning=function(w){NULL})
      if (is.null(coxfit)) {
         P.cph.int[j] <- NA
      } else {
         P.cph.int[j] <- summary(coxfit)$coefficients[3,"Pr(>|z|)",drop=T]
      }
   }
   P.cph.bh.int <- P.cph.int
   nna <- !is.na(P.cph.int)
   P.cph.bh.int[nna] <- p.adjust(p=P.cph.int[nna], method="BH")
   P.cph.bh.int[!nna] <- NA
   
   return(list("raw"=P.cph.int, "fdr"=P.cph.bh.int))
}
#==========================================================================================#



#==========================================================================================#
# Function to display the log file NEWS of updates of the IRSF package.
#==========================================================================================#
IRSF.news <- function(...) {
   
   newsfile <- file.path(system.file(package="IRSF"), "NEWS")
   file.show(newsfile)
   
}
#==========================================================================================#
