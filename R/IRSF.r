#==========================================================================================#
# Ranking of individual and noise variables
# by univariate minimal depth of a maximal subtree (MDMS)
# Dataset X assumes that:
# - all variables are in columns
# - the observed time and censoring variables are in the first two columns
# - each variable has a unique name, excluding the word "noise"
#==========================================================================================#

rsf.rank <- function(X,
                     ntree,
                     method,
                     splitrule="logrank",
                     importance="random",
                     B,
                     ci=90,
                     parallel=FALSE,
                     conf=NULL,
                     verbose=TRUE,
                     seed=NULL) {

   if (!parallel) {
      if (is.null(seed)) {
         digits <- getOption("digits")
         seed <- runif(n=B, min=1, max=2) * 10^(digits-2)
      } else {
         seed <- (0:(B-1)) + seed
      }
      rsf.obj <- rsf.rank.signif(X=X,
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
                              varlist=c("rsf.rank.signif"),
                              envir=.GlobalEnv)
      rsf.cl <- parallel::clusterCall(cl=cl, fun=rsf.rank.signif,
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
# Ranking of interactions between pairs of individual and noise variables by
# bivariate interaction minimal depth of a maximal subtree (IMDMS)
# Dataset X assumes that:
# - all variables are in columns
# - the observed time and censoring variables are in the first two columns
# - each variable has a unique name, excluding the word "noise"
#==========================================================================================#

rsf.int <- function(X,
                    ntree,
                    method,
                    splitrule="logrank",
                    importance="random",
                    B,
                    ci=90,
                    parallel=FALSE,
                    conf=NULL,
                    verbose=TRUE,
                    seed=NULL) {

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
      if (method == "mdms") {
         rsf.obs.boot <- array(data=NA, dim=c(dim(rsf.cl[[1]]$boot.obs)[1], dim(rsf.cl[[1]]$boot.obs)[2], 0))
         rsf.noise.boot <- array(data=NA, dim=c(dim(rsf.cl[[1]]$boot.noise)[1], dim(rsf.cl[[1]]$boot.noise)[2], 0))
         for (b in 1:cpus) {
            rsf.obs.boot <-  abind(rsf.obs.boot, rsf.cl[[b]]$boot.obs)
            rsf.noise.boot <-  abind(rsf.noise.boot, rsf.cl[[b]]$boot.noise)
         }
      } else {
         stop("Unmatched method \n")
      }
   }
   theta <- (100-ci)/200
   if (method == "mdms") {
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
# Subroutine for ranking of individual and noise variables by
# univariate minimal depth of a maximal subtree (MDMS)
# Dataset X assumes that:
# - all variables are in columns
# - the observed time and censoring variables are in the first two columns
# - each variable has a unique name, excluding the word "noise"
#==========================================================================================#

rsf.rank.signif <- function(X,
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
   rank.obs <- matrix(data=NA, nrow=B, ncol=p)
   rank.noise <- matrix(data=NA, nrow=B, ncol=p)
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

      rsf.obs.bo <- rfsrc(formula=Surv(time=time, event=event, type="right") ~ .,
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

      rsf.noise.bo <- rfsrc(formula=Surv(time=time, event=event, type="right") ~ .,
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

      if (method == "mdms") {

         max.obj <- max.subtree(object=rsf.obs.bo,
                                max.order=1,
                                sub.order=TRUE,
                                conservative=FALSE)
         rank.obs[b,] <- diag(max.obj$sub.order)

         max.obj <- max.subtree(object=rsf.noise.bo,
                                max.order=1,
                                sub.order=TRUE,
                                conservative=FALSE)
         rank.noise[b,] <- diag(max.obj$sub.order)

      } else {
         stop("Unmatched method \n")
      }
   }
   colnames(rank.obs) <- rsf.obs.bo$xvar.names
   colnames(rank.noise) <- rsf.noise.bo$xvar.names
   return(list("boot.obs"=rank.obs, "boot.noise"=rank.noise))

}
#==========================================================================================#




#==========================================================================================#
# Subroutine for ranking of interactions between pairs of individual and noise variables by
# bivariate interaction minimal depth of a maximal subtree (IMDMS)
# Dataset X assumes that:
# - all variables are in columns
# - the observed time and censoring variables are in the first two columns
# - each variable has a unique name, excluding the word "noise"
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
   if (method == "mdms") {
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

         rsf.obs.bo <- rfsrc(formula=Surv(time=time, event=event, type="right") ~ .,
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

         rsf.noise.bo <- rfsrc(formula=Surv(time=time, event=event, type="right") ~ .,
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

         int.obs.bo <- find.interaction(object=rsf.obs.bo,
                                        importance=importance,
                                        method="maxsubtree",
                                        sorted=TRUE,
                                        nrep=1,
                                        na.action="na.omit",
                                        seed=seed[b],
                                        verbose=TRUE)

         int.noise.bo <- find.interaction(object=rsf.noise.bo,
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
   } else {
      stop("Unmatched method \n")
   }

}
#==========================================================================================#




#==========================================================================================#
# Computes p-values of significance of regression coefficients of pairwise interaction effects in a Cox-PH model
#==========================================================================================#

cph.int <- function (X, int.term) {

   p <- dim(X)[2] - 2
   int.term.list <- strsplit(x=int.term, split=":")
   fmla.int <- as.formula(paste("Surv(time=", colnames(X)[1], ", event=", colnames(X)[2], ", type=\"right\") ~ . + (.)^2", sep=""))

   P.cph.int <- numeric(p*(p-1)/2)
   names(P.cph.int) <- int.term
   for (j in 1:(p*(p-1)/2)) {
      Z <- X[,c(colnames(X)[1], colnames(X)[2], int.term.list[[j]])]
      coxfit <- tryCatch({coxph(fmla.int, data=Z, model=T, x=T, y=T)}, error=function(w){NULL}, warning=function(w){NULL})
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
# Computes p-values of significance of regression coefficients of main effects in a Cox-PH model
#==========================================================================================#

cph.main <- function (X, main.term) {

   fmla.main <- as.formula(paste("Surv(time=", colnames(X)[1], ", event=", colnames(X)[2], ", type=\"right\") ~ .", sep=""))

   P.cph.main <- numeric(p)
   names(P.cph.main) <- main.term
   for (j in 1:p) {
      Z <- X[,c(colnames(X)[1], colnames(X)[2], main.term[j])]
      coxfit <- tryCatch({coxph(fmla.main, data=Z, model=T, x=T, y=T)}, error=function(w){NULL}, warning=function(w){NULL})
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
# Plot rules of RSF statistics
#==========================================================================================#

plot.rule <- function(x, rule="ci", col, lty, lwd, len=NULL, ...) {

   if (rule == "ci") {

      segments(x0=x[,"obs.LBCI"], y0=x[,"noise.LBCI"],
               x1=x[,"obs.UBCI"], y1=x[,"noise.LBCI"], col=col, lty=lty, lwd=lwd, ...)
      segments(x0=x[,"obs.LBCI"], y0=x[,"noise.UBCI"],
               x1=x[,"obs.UBCI"], y1=x[,"noise.UBCI"], col=col, lty=lty, lwd=lwd, ...)
      segments(x0=x[,"obs.LBCI"], y0=x[,"noise.LBCI"],
               x1=x[,"obs.LBCI"], y1=x[,"noise.UBCI"], col=col, lty=lty, lwd=lwd, ...)
      segments(x0=x[,"obs.UBCI"], y0=x[,"noise.LBCI"],
               x1=x[,"obs.UBCI"], y1=x[,"noise.UBCI"], col=col, lty=lty, lwd=lwd, ...)

   } else if (rule == "1se") {

      arrows(x0=x[,"obs.mean"],
             x1=x[,"obs.mean"]-x[,"obs.se"],
             y0=x[,"noise.mean"],
             y1=x[,"noise.mean"],
             length=len, angle=90, code=2, col=col, lty=lty, lwd=lwd)
      arrows(x0=x[,"obs.mean"],
             x1=x[,"obs.mean"]+x[,"obs.se"],
             y0=x[,"noise.mean"],
             y1=x[,"noise.mean"],
             length=len, angle=90, code=2, col=col, lty=lty, lwd=lwd)

   } else {

      warning("Not a valid rule\n")

   }
}
#==========================================================================================#




#==========================================================================================#
# Overall ensemble Nelson-Aalen estimator
#==========================================================================================#

aalen.estim <- function(obj,
                        subset=NULL) {

   get.event.info <- function(obj, subset) {
      if (grepl("surv", obj$family)) {
         if (!is.null(obj$yvar)) {
            if (is.null(subset)) {
               subset <- (1:nrow(cbind(obj$yvar)))
            }
            r.dim <- 2
            time <- obj$yvar[subset, 1]
            cens <- obj$yvar[subset, 2]
            event <- na.omit(cens)[na.omit(cens) > 0]
            event.type <- unique(event)
         } else {
            r.dim <- 0
            event <- event.type <- cens <- cens <- time <- NULL
         }
         time.interest <- obj$time.interest
      } else {
         r.dim <- 1
         event <- event.type <- cens <- time.interest <- cens <- time <- NULL
      }
      return(list(event=event, event.type=event.type, cens=cens, time.interest=time.interest, time=time, r.dim=r.dim))
   }

   event.info <- get.event.info(obj=obj, subset=subset)
   km.obj <- matrix(unlist(mclapply(X=1:length(event.info$time.interest),
                                    FUN=function(j) {
                                       c(sum(event.info$time >= event.info$time.interest[j], na.rm=TRUE),
                                         sum(event.info$time[event.info$cens != 0] == event.info$time.interest[j], na.rm=TRUE))
                                    }
   )
   ), ncol=2, byrow=TRUE)
   Y <- km.obj[, 1]
   d <- km.obj[, 2]
   r <- d/(Y + 1 * (Y == 0))
   return(exp(-cumsum(r)))

}
#==========================================================================================#




#==========================================================================================#
# Overall ensemble averaged event-free estimator
#==========================================================================================#

avgensb.estim <- function(obj,
                          subset=NULL) {

   if (!is.null(obj$yvar)) {
      if (is.null(subset)) {
         subset <- (1:nrow(cbind(obj$yvar)))
      }
      if  (!is.null(obj$survival.oob)) {
         return(apply(obj$survival.oob[subset,, drop=FALSE], 2, mean, na.rm=TRUE))
      } else {
         return(apply(obj$survival[subset,, drop=FALSE], 2, mean, na.rm=TRUE))
      }
   } else {
      subset <- 1:obj$n
      if  (!is.null(obj$survival.oob)) {
         return(apply(obj$survival.oob[subset,, drop=FALSE], 2, mean, na.rm=TRUE))
      } else {
         return(apply(obj$survival[subset,, drop=FALSE], 2, mean, na.rm=TRUE))
      }
   }

}
#==========================================================================================#
