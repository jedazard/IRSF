###################################################################################################################################
# Correlation of Genomic Variants to Time-To-Event
###################################################################################################################################

#==========================================================================================
# Loading IRSF library for simulations, validations, and publication of the manuscript figures
#==========================================================================================
library("IRSF")

#==========================================================================================
# Loading additional libraries
#==========================================================================================
library("ggRandomForests")
library("NADA")
library("MASS")
library("RColorBrewer")



##########################################################################################################################################
# For constants, parameters and functions
##########################################################################################################################################

#=========================================================================================#
# Cluster configuration
#=========================================================================================#
if (require("parallel")) {
  print("'parallel' is attached correctly \n")
} else {
  stop("'parallel' must be attached first \n")
}

if (.Platform$OS.type == "windows") {
  # Windows OS
  cpus <- detectCores(logical = TRUE)
  conf <- list("spec"=rep("localhost", cpus),
               "type"="SOCKET",
               "homo"=TRUE,
               "verbose"=TRUE,
               "outfile"=file.path(getwd(), "output_SOCK.txt", fsep=.Platform$file.sep))
} else if (.Platform$OS.type == "unix") {
  # Linux or Mac OS
  argv <- commandArgs(trailingOnly=TRUE)  
  if (is.empty(argv)) {
    # Mac OS : No argument "argv" in this case  
    cpus <- detectCores(logical = TRUE)
    conf <- list("spec"=rep("localhost", cpus),
                 "type"="SOCKET",
                 "homo"=TRUE,
                 "verbose"=TRUE,
                 "outfile"=file.path(getwd(), "output_SOCK.txt", fsep=.Platform$file.sep))
  } else {
    # Linux OS : Retrieving arguments "type" and "cpus" from the SLURM script
    type <- as.character(argv[1])
    cpus <- as.numeric(Sys.getenv("SLURM_NTASKS"))
    if (type == "SOCKET") {
      conf <- list("spec"=rep("localhost", cpus),
                   "type"="SOCKET",
                   "homo"=TRUE,
                   "verbose"=TRUE,
                   "outfile"=file.path(getwd(), "output_SOCK.txt", fsep=.Platform$file.sep))
    } else if (type == "MPI") {
      if (require("Rmpi")) {
          print("Rmpi is loaded correctly \n")
      } else {
          stop("Rmpi must be installed first to use MPI\n")
      }
      conf <- list("spec"=cpus,
                   "type"="MPI",
                   "homo"=TRUE,
                   "verbose"=TRUE,
                   "outfile"=file.path(getwd(), "output_MPI.txt", fsep=.Platform$file.sep))
    } else {
      stop("Unrecognized cluster type: you must specify a \"SOCKET\" or \"MPI\" cluster type\n")
    }
  }
} else {
  stop("Unrecognized platform: you must specify a \"windows\" or \"unix\" platform type\n")
}

cat("Cluster configuration:\n")
print(conf)


##########################################################################################################################################
# Real Dataset Study
# Default parameters were used for all main functions.
##########################################################################################################################################

#=========================================================================================#
# Constants
#=========================================================================================#
ntree <- 1000

#============================================================================================
# Read in the data
#============================================================================================
data("MACS", package="IRSF")

formula.X <- as.formula(Surv(time=TTX, event=EventX, type="right") ~ 1  + Race + Group3 + DEFB.CNV3 + CCR2.SNP + CCR5.SNP2 + CCR5.ORF + CXCL12.SNP2)
formula.A <- as.formula(Surv(time=TTA, event=EventA, type="right") ~ 1  + Race + Group3 + DEFB.CNV3 + CCR2.SNP + CCR5.SNP2 + CCR5.ORF + CXCL12.SNP2)

X <- MACS[,c("TTX","EventX","Race","Group3","DEFB.CNV3","CCR2.SNP","CCR5.SNP2","CCR5.ORF","CXCL12.SNP2")]
A <- MACS[,c("TTA","EventA","Race","Group3","DEFB.CNV3","CCR2.SNP","CCR5.SNP2","CCR5.ORF","CXCL12.SNP2")]

n <- nrow(X)
p <- ncol(X) - 2

#==========================================================================================
# Global error rate from Random Survival Forests
# The prediction error rate is based on the predicted value of mortality, which indirectly arises from the CHF
#==========================================================================================

#=====================
# X4-Emergence outcome
#=====================
allX.rfsrc <- rfsrc(formula=formula.X,
                    data=X,
                    ntree=ntree,
                    bootstrap="by.root",
                    mtry=p,
                    nodesize=3,
                    splitrule="logrank",
                    nsplit=0,
                    importance="random",
                    na.action="na.omit",
                    proximity=TRUE,
                    samptype="swr",
                    forest=TRUE,
                    var.used="all.trees",
                    split.depth="all.trees",
                    membership=TRUE,
                    statistics=TRUE,
                    tree.err=TRUE,
                    seed=12345678)

allX.rfsrc
                    Sample size: 50
                    Number of deaths: 32
                     Number of trees: 1000
          Minimum terminal node size: 3
       Average no. of terminal nodes: 14.202
No. of variables tried at each split: 7
              Total no. of variables: 7
                            Analysis: RSF
                              Family: surv
                      Splitting rule: logrank
                          Error rate: 50.53%

summary(allX.rfsrc)
allX.rfsrc$err.rate         # ensemble RSF cumulative OOB error rate, i.e. based on OOB data.
allX.rfsrc$importance       # ensemble Variable importance (VIMP) for each x-variable
allX.rfsrc$predicted.oob    # ensemble OOB predicted eventuality (X4-Emergence) value
allX.rfsrc$survival.oob     # ensemble OOB event-free (survival) function
allX.rfsrc$chf.oob          # ensemble OOB cumulative hazard function (CHF)
allX.rfsrc$time.interest    # ordered unique event times

#========================
# AIDS-Diagnosis outcome
#========================
allA.rfsrc <- rfsrc(formula=formula.A,
                    data=A,
                    ntree=ntree,
                    bootstrap="by.root",
                    mtry=p,
                    nodesize=3,
                    splitrule="logrank",
                    nsplit=0,
                    importance="random",
                    na.action="na.omit",
                    proximity=TRUE,
                    samptype="swr",
                    forest=TRUE,
                    var.used="all.trees",
                    split.depth="all.trees",
                    membership=TRUE,
                    statistics=TRUE,
                    tree.err=TRUE,
                    seed=12345678)

allA.rfsrc
                    Sample size: 50
                    Number of deaths: 49
                     Number of trees: 1000
          Minimum terminal node size: 3
       Average no. of terminal nodes: 14.453
No. of variables tried at each split: 7
              Total no. of variables: 7
                            Analysis: RSF
                              Family: surv
                      Splitting rule: logrank
                          Error rate: 32.32%

summary(allA.rfsrc)
allA.rfsrc$err.rate         # ensemble RSF cumulative OOB error rate, i.e. based on OOB data.
allA.rfsrc$importance       # ensemble Variable importance (VIMP) for each x-variable
allA.rfsrc$predicted.oob    # ensemble OOB predicted eventuality (AIDS-Diagnosis) value
allA.rfsrc$survival.oob     # ensemble OOB event-free (survival) function
allA.rfsrc$chf.oob          # ensemble OOB cumulative hazard function (CHF)
allA.rfsrc$time.interest    # ordered unique event times


#==========================================================================================
# Prediction estimates
#==========================================================================================
set.seed(12345678)
train <- sample(1:n, round(n * 3/5))
test <- setdiff(1:n, train)

#=====================
# primary call
#=====================
allX.train.rfsrc <- rfsrc(formula=formula.X,
                          data=X[train, ],
                          ntree=ntree,
                          bootstrap="by.root",
                          mtry=p,
                          nodesize=3,
                          splitrule="logrank",
                          nsplit=0,
                          importance="random",
                          na.action="na.omit",
                          proximity=TRUE,
                          samptype="swr",
                          forest=TRUE,
                          var.used="all.trees",
                          split.depth="all.trees",
                          membership=TRUE,
                          statistics=TRUE,
                          tree.err=TRUE,
                          seed=12345678)

allA.train.rfsrc <- rfsrc(formula=formula.A,
                          data=A[train, ],
                          ntree=ntree,
                          bootstrap="by.root",
                          mtry=p,
                          nodesize=3,
                          splitrule="logrank",
                          nsplit=0,
                          importance="random",
                          na.action="na.omit",
                          proximity=TRUE,
                          samptype="swr",
                          forest=TRUE,
                          var.used="all.trees",
                          split.depth="all.trees",
                          membership=TRUE,
                          statistics=TRUE,
                          tree.err=TRUE,
                          seed=12345678)

#=====================
# Non-Standard prediction estimates:
# overlays the 'test' data on the train-grown forest
#=====================
allX.pred <- predict(object=allX.train.rfsrc,
                     newdata=X[test, ],
                     importance="random",
                     na.action="na.omit",
                     outcome="test",
                     proximity=TRUE,
                     var.used="all.trees",
                     split.depth="all.trees",
                     membership=TRUE,
                     statistics=TRUE,
                     do.trace=TRUE,
                     seed=12345678)

allX.pred

allA.pred <- predict(object=allA.train.rfsrc,
                     newdata=A[test, ],
                     importance="random",
                     na.action="na.omit",
                     outcome="test",
                     proximity=TRUE,
                     var.used="all.trees",
                     split.depth="all.trees",
                     membership=TRUE,
                     statistics=TRUE,
                     do.trace=TRUE,
                     seed=12345678)

allA.pred


#==========================================================================================
# Plot of global prediction performance and Variable Importance from a RF-SRC analysis
# Cumulative OOB Error Rate is between 0 and 1, and measures how well the predictor correctly
# ranks (classifies) random individuals in terms of event probability.
# A value of 0.5 is no better than random guessing.
# A value of 0 is perfect.
#==========================================================================================
vimpX.obj <- vimp(object=allX.rfsrc,
                  importance="random",
                  joint=FALSE,
                  seed=12345678)

vimpA.obj <- vimp(object=allA.rfsrc,
                  importance="random",
                  joint=FALSE,
                  seed=12345678)

ooberX <- sort(vimpX.obj$importance, decreasing=T)
ooberA <- sort(vimpA.obj$importance, decreasing=T)


#==========================================================================================
# Ranking of individual and noise variables by
# univariate minimal depth of a maximal subtree (MDMS)
#==========================================================================================
#====
# TTX
#====
allX.main.mdms <- rsf.main(X=MACS[,c("TTX","EventX","Race","Group3","DEFB.CNV3","CCR2.SNP","CCR5.SNP2","CCR5.ORF","CXCL12.SNP2")],
                           ntree=ntree, method="mdms", splitrule="logrank", importance="random", B=1000, ci=80,
                           parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)
#====
# TTA
#====
allA.main.mdms <- rsf.main(X=MACS[,c("TTA","EventA","Race","Group3","DEFB.CNV3","CCR2.SNP","CCR5.SNP2","CCR5.ORF","CXCL12.SNP2")],
                           ntree=ntree, method="mdms", splitrule="logrank", importance="random", B=1000, ci=80,
                           parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

#========================================================================================#
# Fits a Proportional Hazards Time-To-Event Regression Model saturated with first order terms.
# Computes p-values of significance of regression coefficients of main effects in a Cox-PH model
#==========================================================================================#
allX.main.cph <- cph.main(X=MACS[,c("TTX","EventX","Race","Group3","DEFB.CNV3","CCR2.SNP","CCR5.SNP2","CCR5.ORF","CXCL12.SNP2")], main.term=rownames(allX.main.mdms))
allA.main.cph <- cph.main(X=MACS[,c("TTA","EventA","Race","Group3","DEFB.CNV3","CCR2.SNP","CCR5.SNP2","CCR5.ORF","CXCL12.SNP2")], main.term=rownames(allA.main.mdms))

#==========================================================================================
# Ranking of interactions between pairs of individual and noise variables by
# bivariate interaction minimal depth of a maximal subtree (IMDMS)
# Entries [i][j] indicate the normalized minimal depth of a variable [j] w.r.t. the maximal subtree for variable [i]
# (normalized w.r.t. the size of [i]'s maximal subtree).
#==========================================================================================
#====
# TTX
#====
allX.int.mdms <- rsf.int(X=MACS[,c("TTX","EventX","Race","Group3","DEFB.CNV3","CCR2.SNP","CCR5.SNP2","CCR5.ORF","CXCL12.SNP2")],
                         ntree=ntree, method="imdms", splitrule="logrank", importance="random", B=1000, ci=80,
                         parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

#====
# TTA
#====
allA.int.mdms <- rsf.int(X=MACS[,c("TTA","EventA","Race","Group3","DEFB.CNV3","CCR2.SNP","CCR5.SNP2","CCR5.ORF","CXCL12.SNP2")],
                         ntree=ntree, method="imdms", splitrule="logrank", importance="random", B=1000, ci=80,
                         parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

#==========================================================================================
# Fits a Proportional Hazards Time-To-Event Regression Model saturated with first and second order terms
# Computes p-values of significance of regression coefficients of pairwise interaction effects in a Cox-PH model
#==========================================================================================
allX.int.cph <- cph.int(X=MACS[,c("TTX","EventX","Race","Group3","DEFB.CNV3","CCR2.SNP","CCR5.SNP2","CCR5.ORF","CXCL12.SNP2")], int.term=rownames(allX.int.mdms))
allA.int.cph <- cph.int(X=MACS[,c("TTA","EventA","Race","Group3","DEFB.CNV3","CCR2.SNP","CCR5.SNP2","CCR5.ORF","CXCL12.SNP2")], int.term=rownames(allA.int.mdms))



