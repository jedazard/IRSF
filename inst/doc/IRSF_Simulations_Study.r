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
# Simulations Study
# Default parameters were used for all main functions.
##########################################################################################################################################

#==========================================================================================#
# Constants
#==========================================================================================#
ntree <- 1000

#==========================================================================================#
# Continuous case:
# All variables xj, j in {1,...,p}, are iid from a multivariate uniform distribution
# with parmeters  a=1, b=5, i.e. on [1, 5].
# rho = 0.50
#==========================================================================================#
seed <- 1234567
set.seed(seed)

n <- 200
p <- 5
x <- matrix(data=runif(n=n*p, min=1, max=5),
            nrow=n, ncol=p, byrow=FALSE,
            dimnames=list(1:n, paste("X", 1:p, sep="")))

#===================================================#
# negative control
#===================================================#
#neg 0
# ----
beta <- rep(0,p)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.52)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg0 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.95)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg0 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#neg 1
# X1 + X2 + X3 + X4 + X5
beta <- rep(1,p)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.33)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg1 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=31)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg1 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#neg 2
# X1 + X5
beta <- c(1,0,0,0,1)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.50)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg2 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=13.8)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg2 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#===================================================#
# positive controls
#===================================================#
#pos 1
# X1X2
beta <- c(rep(0,p), 1)
covar <- cbind(x, "X1X2"=x[,1]*x[,2])
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=3.33)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos1 <- data.frame(stime, status, x)
summary(status)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=20.80)                            # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos1 <- data.frame(stime, status, x)
summary(status)

#pos 2
# X1 + X2 + X1X2
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(x, "X1X2"=x[,1]*x[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=3.9)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=33.8)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2 <- data.frame(stime, status, x)
summary(status)

#pos 3
# X1 + X2 + X3 + X4 + X5 + X1X2
beta <- c(rep(1,p), 1)
covar <- cbind(x, "X1X2"=x[,1]*x[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=3.8)                                # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos3 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=52)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos3 <- data.frame(stime, status, x)
summary(status)

#pos 4
# X1 + X2 + exp(-(X1X2))
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(x, "exp(-(X1X2))"=exp(-(x[,1]*x[,2])))
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.4)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos4 <- data.frame(stime, status, x)
summary(status)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=14.2)                              # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos4 <- data.frame(stime, status, x)
summary(status)

#pos 5
# X1 + X2 + sqrt(X1X2)
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(x, "sqrt(X1X2)"=sqrt(x[,1]*x[,2]))
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.6)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos5 <- data.frame(stime, status, x)
summary(status)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=20.22)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos5 <- data.frame(stime, status, x)
summary(status)

#pos 6
# X1 + X2 + log(X1X2)
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(x, "log(X1X2)"=log(x[,1]*x[,2]))
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.5)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos6 <- data.frame(stime, status, x)
summary(status)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=18.5)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos6 <- data.frame(stime, status, x)
summary(status)

#pos 7
# X1 + X2 + sin(pi.X1X2)
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(x, "sin(piX1X2)"=sin(pi*x[,1]*x[,2]))
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.48)                                # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos7 <- data.frame(stime, status, x)
summary(status)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=13.90)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos7 <- data.frame(stime, status, x)
summary(status)


#==========================================================================================#
# Effect of censoring
# Continuous Case
# regression model #pos 2
# X1 + X2 + X1X2
#==========================================================================================#
#pos 2.1
# rho = 0.30
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(x, "X1X2"=x[,1]*x[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=0.08)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.1 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=22)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.1 <- data.frame(stime, status, x)
summary(status)

#pos 2.2
# rho = 0.50
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(x, "X1X2"=x[,1]*x[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=3.9)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.2 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=33.8)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.2 <- data.frame(stime, status, x)
summary(status)

#pos 2.3
# rho = 0.70
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(x, "X1X2"=x[,1]*x[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=89)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.3 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=57)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.3 <- data.frame(stime, status, x)
summary(status)


#==========================================================================================#
# Effect of variable
# Continuous Case
# rho = 0.50
#==========================================================================================#
#neg 2.1
# X1 + X2
beta <- c(1,1,0,0,0)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.35)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg2.1 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=14.2)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg2.1 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#neg 2.2
# X1 + X3
beta <- c(1,0,1,0,0)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.4)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg2.2 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=13.85)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg2.2 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#neg 2.3
# X1 + X4
beta <- c(1,0,0,1,0)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.39)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg2.3 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=14.2)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg2.3 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#neg 2.4
# X1 + X5
beta <- c(1,0,0,0,1)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.50)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg2.4 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=13.8)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg2.4 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#neg 2.5
# X2 + X3
beta <- c(0,1,1,0,0)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.7)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg2.5 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=13.5)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg2.5 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#neg 2.6
# X2 + X4
beta <- c(0,1,0,1,0)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.72)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg2.6 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=14.4)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg2.6 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#neg 2.7
# X2 + X5
beta <- c(0,1,0,0,1)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.9)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg2.7 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=14.3)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg2.7 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#neg 2.8
# X3 + X4
beta <- c(0,0,1,1,0)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.60)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg2.8 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=13.5)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg2.8 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#neg 2.9
# X3 + X5
beta <- c(0,0,1,0,1)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.257)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg2.9 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=13.8)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg2.9 <- data.frame(stime, status, x)
summary(status)
summary(stime)

#neg 2.10
# X4 + X5
beta <- c(0,0,0,1,1)
covar <- x
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.72)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.neg2.10 <- data.frame(stime, status, x)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=13.8)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.neg2.10 <- data.frame(stime, status, x)
summary(status)
summary(stime)



#pos 2.4
# X1 + X2 + X1X2
beta <- c(1, 1, 0, 0, 0, 1)
covar <- cbind(x, "X1X2"=x[,1]*x[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=3.9)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.4 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=33.8)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.4 <- data.frame(stime, status, x)
summary(status)

#pos 2.5
# X1 + X3 + X1X3
beta <- c(1, 0, 1, 0, 0, 1)
covar <- cbind(x, "X1X3"=x[,1]*x[,3])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=2.05)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.5 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=34.5)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.5 <- data.frame(stime, status, x)
summary(status)

#pos 2.6
# X1 + X4 + X1X4
beta <- c(1, 0, 0, 1, 0, 1)
covar <- cbind(x, "X1X4"=x[,1]*x[,4])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=2.67)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.6 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=32.85)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.6 <- data.frame(stime, status, x)
summary(status)

#pos 2.7
# X1 + X5 + X1X5
beta <- c(1, 0, 0, 0, 1, 1)
covar <- cbind(x, "X1X5"=x[,1]*x[,5])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=2.1)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.7 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=32.5)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.7 <- data.frame(stime, status, x)
summary(status)

#pos 2.8
# X2 + X3 + X2X3
beta <- c(0, 1, 1, 0, 0, 1)
covar <- cbind(x, "X2X3"=x[,2]*x[,3])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=3.15)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.8 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=36.07)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.8 <- data.frame(stime, status, x)
summary(status)

#pos 2.9
# X2 + X4 + X2X4
beta <- c(0, 1, 0, 1, 0, 1)
covar <- cbind(x, "X2X4"=x[,2]*x[,4])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=2.90)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.9 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=33.95)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.9 <- data.frame(stime, status, x)
summary(status)

#pos 2.10
# X2 + X5 + X2X5
beta <- c(0, 1, 0, 0, 1, 1)
covar <- cbind(x, "X2X5"=x[,2]*x[,5])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=3.7)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.10 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=37)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.10 <- data.frame(stime, status, x)
summary(status)

#pos 2.11
# X3 + X4 + X3X4
beta <- c(0, 0, 1, 1, 0, 1)
covar <- cbind(x, "X3X4"=x[,3]*x[,4])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=4.8)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.11 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=32.5)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.11 <- data.frame(stime, status, x)
summary(status)

#pos 2.12
# X3 + X5 + X3X5
beta <- c(0, 0, 1, 0, 1, 1)
covar <- cbind(x, "X3X5"=x[,3]*x[,5])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=2.80)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.12 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=35.192)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.12 <- data.frame(stime, status, x)
summary(status)

#pos 2.13
# X4 + X5 + X4X5
beta <- c(0, 0, 0, 1, 1, 1)
covar <- cbind(x, "X4X5"=x[,4]*x[,5])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=5.00)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cont.pos2.13 <- data.frame(stime, status, x)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=33.2)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cont.pos2.13 <- data.frame(stime, status, x)
summary(status)


rm(x, covar, beta, eta0, eta, lambda0, lambda, tt, tc, stime, status)


#==========================================================================================#
# Discrete case (5 categories):
# All variables xj, j in {1,...,p}, are iid from a multivariate binomial distribution
# with equiprobability p=1/5.
# rho = 0.50
#==========================================================================================#
seed <- 1234567
set.seed(seed)

n <- 200
p <- 5
z <- matrix(data=sample(c(1,2,3,4,5), size=n*p, replace=TRUE, prob=c(1/5,1/5,1/5,1/5,1/5)),
            nrow=n, ncol=p, byrow=FALSE,
            dimnames=list(1:n, paste("X", 1:p, sep="")))

# z.scale <- scale(x=z, center=FALSE, scale=apply(z, 2, sd, na.rm = TRUE))

#===================================================#
# negative controls
#===================================================#
#neg 0
# ----
beta <- rep(0,p)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.52)                                # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg0 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.95)                                # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg0 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#neg 1
# X1 + X2 + X3 + X4 + X5
beta <- rep(1,p)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.67)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg1 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=35.1)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg1 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#neg 2
# X1 + X5
beta <- c(1,0,0,0,1)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.41)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg2 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=16.503)                           # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg2 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#===================================================#
# positive controls
#===================================================#
#pos 1
# X1X2
beta <- c(rep(0,p), 1)
covar <- cbind(z, "X1X2"=z[,1]*z[,2])
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=3.5)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos1 <- data.frame(stime, status, z)
summary(status)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=26.3)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos1 <- data.frame(stime, status, z)
summary(status)

#pos 2
# X1 + X2 + X1X2
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(z, "X1X2"=z[,1]*z[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=5.2)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=33)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2 <- data.frame(stime, status, z)
summary(status)

#pos 3
# X1 + X2 + X3 + X4 + X5 + X1X2
beta <- c(rep(1,p), 1)
covar <- cbind(z, "X1X2"=z[,1]*z[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=5.5)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos3 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=52)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos3 <- data.frame(stime, status, z)
summary(status)

#pos 4
# X1 + X2 + exp(-(X1X2))
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(z, "exp(-(X1X2))"=exp(-(z[,1]*z[,2])))
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.72)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos4 <- data.frame(stime, status, z)
summary(status)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=16.2)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos4 <- data.frame(stime, status, z)
summary(status)

#pos 5
# X1 + X2 + sqrt(X1X2)
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(z, "sqrt(X1X2)"=sqrt(z[,1]*z[,2]))
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.55)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos5 <- data.frame(stime, status, z)
summary(status)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=23.122)                           # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos5 <- data.frame(stime, status, z)
summary(status)

#pos 6
# X1 + X2 + log(X1X2)
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(z, "log(X1X2)"=log(z[,1]*z[,2]))
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.365)                            # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos6 <- data.frame(stime, status, z)
summary(status)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=20.950)                           # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos6 <- data.frame(stime, status, z)
summary(status)

#pos 7
# X1 + X2 + sin(pi.X1X2)
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(z, "sin(piX1X2)"=sin(pi*z[,1]*z[,2]))
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.7)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos7 <- data.frame(stime, status, z)
summary(status)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=16.1)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos7 <- data.frame(stime, status, z)
summary(status)


#==========================================================================================#
# Effect of censoring
# Discrete Case
# regression model #pos 2
# X1 + X2 + X1X2
#==========================================================================================#
#pos 2.1
# rho = 0.30
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(z, "X1X2"=z[,1]*z[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=0.035)                            # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.1 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=21)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.1 <- data.frame(stime, status, z)
summary(status)

#pos 2.2
# rho = 0.50
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(z, "X1X2"=z[,1]*z[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=5.5)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.2 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=33)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.2 <- data.frame(stime, status, z)
summary(status)

#pos 2.3
# rho = 0.70
beta <- c(rep(1,2), rep(0,p-2), 1)
covar <- cbind(z, "X1X2"=z[,1]*z[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=195)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.3 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=57.5)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.3 <- data.frame(stime, status, z)
summary(status)


#==========================================================================================#
# Effect of variable
# Discrete Case
# rho = 0.50
#==========================================================================================#
#neg 2.1
# X1 + X4
beta <- c(1,1,0,0,0)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.65)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg2.1 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=16.1)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg2.1 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#neg 2.2
# X1 + X4
beta <- c(1,0,1,0,0)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.25)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg2.2 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=16.1)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg2.2 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#neg 2.3
# X1 + X4
beta <- c(1,0,0,1,0)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.46)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg2.3 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=16.65)                            # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg2.3 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#neg 2.4
# X1 + X5
beta <- c(1,0,0,0,1)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.41)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg2.4 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=16.503)                           # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg2.4 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#neg 2.5
# X2 + X3
beta <- c(0,1,1,0,0)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.70)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg2.5 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=13.6)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg2.5 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#neg 2.6
# X2 + X4
beta <- c(0,1,0,1,0)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=2.1)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg2.6 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=14.0)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg2.6 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#neg 2.7
# X2 + X5
beta <- c(0,1,0,0,1)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.41)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg2.7 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=14.1)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg2.7 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#neg 2.8
# X3 + X4
beta <- c(0,0,1,1,0)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.25)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg2.8 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=14.3)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg2.8 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#neg 2.9
# X3 + X5
beta <- c(0,0,1,0,1)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.7)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg2.9 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=14.5)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg2.9 <- data.frame(stime, status, z)
summary(status)
summary(stime)

#neg 2.10
# X4 + X5
beta <- c(0,0,0,1,1)
covar <- z
eta <- covar %*% beta                                         # regression function

seed <- 1234567
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=1.55)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.neg2.10 <- data.frame(stime, status, z)
summary(status)
summary(stime)

seed <- 1234567
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=14.1)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.neg2.10 <- data.frame(stime, status, z)
summary(status)
summary(stime)



#pos 2.4
# X1 + X2 + X1X2
beta <- c(1, 1, 0, 0, 0, 1)
covar <- cbind(z, "X1X2"=z[,1]*z[,2])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=5.2)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.4 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=33)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.4 <- data.frame(stime, status, z)
summary(status)

#pos 2.5
# X1 + X3 + X1X3
beta <- c(1, 0, 1, 0, 0, 1)
covar <- cbind(z, "X1X3"=z[,1]*z[,3])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=7.0)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.5 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=34.3)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.5 <- data.frame(stime, status, z)
summary(status)

#pos 2.6
# X1 + X4 + X1X4
beta <- c(1, 0, 0, 1, 0, 1)
covar <- cbind(z, "X1X4"=z[,1]*z[,4])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=4.80)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.6 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=37)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.6 <- data.frame(stime, status, z)
summary(status)

#pos 2.7
# X1 + X5 + X1X5
beta <- c(1, 0, 0, 0, 1, 1)
covar <- cbind(z, "X1X5"=z[,1]*z[,5])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=7.5)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.7 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=36)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.7 <- data.frame(stime, status, z)
summary(status)

#pos 2.8
# X2 + X3 + X2X3
beta <- c(0, 1, 1, 0, 0, 1)
covar <- cbind(z, "X2X3"=z[,2]*z[,3])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=8.9)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.8 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=30.3)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.8 <- data.frame(stime, status, z)
summary(status)

#pos 2.9
# X2 + X4 + X2X4
beta <- c(0, 1, 0, 1, 0, 1)
covar <- cbind(z, "X2X4"=z[,2]*z[,4])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=9.55)                             # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.9 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=33.6)                             # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.9 <- data.frame(stime, status, z)
summary(status)

#pos 2.10
# X2 + X5 + X2X5
beta <- c(0, 1, 0, 0, 1, 1)
covar <- cbind(z, "X2X5"=z[,2]*z[,5])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=3.0)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.10 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=35)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.10 <- data.frame(stime, status, z)
summary(status)

#pos 2.11
# X3 + X4 + X3X4
beta <- c(0, 0, 1, 1, 0, 1)
covar <- cbind(z, "X3X4"=z[,3]*z[,4])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=3.8)                              # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.11 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=34.5)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.11 <- data.frame(stime, status, z)
summary(status)

#pos 2.12
# X3 + X5 + X3X5
beta <- c(0, 0, 1, 0, 1, 1)
covar <- cbind(z, "X3X5"=z[,3]*z[,5])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=6.2)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.12 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=35.3)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.12 <- data.frame(stime, status, z)
summary(status)

#pos 2.13
# X4 + X5 + X4X5
beta <- c(0, 0, 0, 1, 1, 1)
covar <- cbind(z, "X4X5"=z[,4]*z[,5])
eta <- covar %*% beta                                         # regression function

seed <- 123456789
set.seed(seed)
lambda0 <- 1
lambda <- lambda0 * exp(eta - mean(eta))                      # hazards function
tt <- rexp(n=n, rate=lambda)                                  # true (uncensored) failure/event times
tc <- runif(n=n, min=0, max=4.8)                               # true (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed failure/event times
status <- 1 * (tt <= tc)                                      # observed failure/event indicator
exp.cat.pos2.13 <- data.frame(stime, status, z)
summary(status)

seed <- 123456789
set.seed(seed)
eta0 <- 1
tt <- as.vector(eta0 + eta)                                   # true latent (uncensored) failure/event times
tc <- runif(n=n, min=0, max=36.5)                               # true latent (censored) failure/event times
stime <- pmin(tt, tc)                                         # observed latent failure/event times
status <- 1 * (tt <= tc)                                      # observed latent failure/event indicator
llv.cat.pos2.13 <- data.frame(stime, status, z)
summary(status)


rm(z, covar, beta, eta0, eta, lambda0, lambda, tt, tc, stime, status)


#==========================================================================================#
# RSF
# Model #1 (#neg0)
# Model #3 (#neg2)
# rho = 0.50
#==========================================================================================#
formula.exp.cont.neg0 <- as.formula(paste("Surv(time=", colnames(exp.cont.neg0)[1], ", event=", colnames(exp.cont.neg0)[2], ", type=\"right\") ~ .", sep=""))
formula.exp.cat.neg0 <- as.formula(paste("Surv(time=", colnames(exp.cat.neg0)[1], ", event=", colnames(exp.cat.neg0)[2], ", type=\"right\") ~ .", sep=""))
formula.llv.cont.neg0 <- as.formula(paste("Surv(time=", colnames(llv.cont.neg0)[1], ", event=", colnames(llv.cont.neg0)[2], ", type=\"right\") ~ .", sep=""))
formula.llv.cat.neg0 <- as.formula(paste("Surv(time=", colnames(llv.cat.neg0)[1], ", event=", colnames(llv.cat.neg0)[2], ", type=\"right\") ~ .", sep=""))

formula.exp.cont.neg2 <- as.formula(paste("Surv(time=", colnames(exp.cont.neg2)[1], ", event=", colnames(exp.cont.neg2)[2], ", type=\"right\") ~ .", sep=""))
formula.exp.cat.neg2 <- as.formula(paste("Surv(time=", colnames(exp.cat.neg2)[1], ", event=", colnames(exp.cat.neg2)[2], ", type=\"right\") ~ .", sep=""))
formula.llv.cont.neg2 <- as.formula(paste("Surv(time=", colnames(llv.cont.neg2)[1], ", event=", colnames(llv.cont.neg2)[2], ", type=\"right\") ~ .", sep=""))
formula.llv.cat.neg2 <- as.formula(paste("Surv(time=", colnames(llv.cat.neg2)[1], ", event=", colnames(llv.cat.neg2)[2], ", type=\"right\") ~ .", sep=""))

rfsrc.exp.cont.neg0 <- rfsrc(formula=formula.exp.cont.neg0,
                             data=exp.cont.neg0,
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

rfsrc.exp.cat.neg0 <- rfsrc(formula=formula.exp.cat.neg0,
                            data=exp.cat.neg0,
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

rfsrc.llv.cont.neg0 <- rfsrc(formula=formula.llv.cont.neg0,
                             data=llv.cont.neg0,
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

rfsrc.llv.cat.neg0 <- rfsrc(formula=formula.llv.cat.neg0,
                            data=llv.cat.neg0,
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


rfsrc.exp.cont.neg2 <- rfsrc(formula=formula.exp.cont.neg2,
                             data=exp.cont.neg2,
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

rfsrc.exp.cat.neg2 <- rfsrc(formula=formula.exp.cat.neg2,
                            data=exp.cat.neg2,
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

rfsrc.llv.cont.neg2 <- rfsrc(formula=formula.llv.cont.neg2,
                             data=llv.cont.neg2,
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

rfsrc.llv.cat.neg2 <- rfsrc(formula=formula.llv.cat.neg2,
                            data=llv.cat.neg2,
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


#==========================================================================================
# Prediction estimates
#==========================================================================================
set.seed(12345678)
train <- sample(1:n, round(n * 4/5))
test <- setdiff(1:n, train)

#=====================
# primary call
#=====================
rfsrc.train.exp.cont.neg0 <- rfsrc(formula=formula.exp.cont.neg0,
                                   data=exp.cont.neg0[train, ],
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

rfsrc.train.exp.cat.neg0 <- rfsrc(formula=formula.exp.cat.neg0,
                                  data=exp.cat.neg0[train, ],
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

rfsrc.train.llv.cont.neg0 <- rfsrc(formula=formula.llv.cont.neg0,
                                   data=llv.cont.neg0[train, ],
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

rfsrc.train.llv.cat.neg0 <- rfsrc(formula=formula.llv.cat.neg0,
                                  data=llv.cat.neg0[train, ],
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

rfsrc.train.exp.cont.neg2 <- rfsrc(formula=formula.exp.cont.neg2,
                                   data=exp.cont.neg2[train, ],
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

rfsrc.train.exp.cat.neg2 <- rfsrc(formula=formula.exp.cat.neg2,
                                  data=exp.cat.neg2[train, ],
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

rfsrc.train.llv.cont.neg2 <- rfsrc(formula=formula.llv.cont.neg2,
                                   data=llv.cont.neg2[train, ],
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

rfsrc.train.llv.cat.neg2 <- rfsrc(formula=formula.llv.cat.neg2,
                                  data=llv.cat.neg2[train, ],
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
rfsrc.pred.exp.cont.neg0 <- predict(object=rfsrc.train.exp.cont.neg0,
                                    newdata=exp.cont.neg0[test, ],
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

rfsrc.pred.exp.cat.neg0 <- predict(object=rfsrc.train.exp.cat.neg0,
                                    newdata=exp.cat.neg0[test, ],
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

rfsrc.pred.llv.cont.neg0 <- predict(object=rfsrc.train.llv.cont.neg0,
                                    newdata=llv.cont.neg0[test, ],
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

rfsrc.pred.llv.cat.neg0 <- predict(object=rfsrc.train.llv.cat.neg0,
                                   newdata=llv.cat.neg0[test, ],
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

rfsrc.pred.exp.cont.neg2 <- predict(object=rfsrc.train.exp.cont.neg2,
                                    newdata=exp.cont.neg2[test, ],
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

rfsrc.pred.exp.cat.neg2 <- predict(object=rfsrc.train.exp.cat.neg2,
                                   newdata=exp.cat.neg2[test, ],
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

rfsrc.pred.llv.cont.neg2 <- predict(object=rfsrc.train.llv.cont.neg2,
                                    newdata=llv.cont.neg2[test, ],
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

rfsrc.pred.llv.cat.neg2 <- predict(object=rfsrc.train.llv.cat.neg2,
                                   newdata=llv.cat.neg2[test, ],
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

#==========================================================================================#
# Ranking of individual and noise variables main effects 
# by univariate minimal depth of a maximal subtree (MDMS)
# Positive and Negative Controls
# rho = 0.50
#==========================================================================================#
main.mdms.exp.cont.neg0 <- rsf.main(X=exp.cont.neg0, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                    B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cat.neg0 <- rsf.main(X=exp.cat.neg0, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                   B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cont.neg0 <- rsf.main(X=llv.cont.neg0, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                    B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cat.neg0 <- rsf.main(X=llv.cat.neg0, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                   B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)


main.mdms.exp.cont.neg2 <- rsf.main(X=exp.cont.neg2, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                    B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cat.neg2 <- rsf.main(X=exp.cat.neg2, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                   B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cont.neg2 <- rsf.main(X=llv.cont.neg2, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                    B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cat.neg2 <- rsf.main(X=llv.cat.neg2, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                   B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)


#==========================================================================================#
# Ranking of individual and noise variables main effects 
# by univariate minimal depth of a maximal subtree (MDMS)
# Effect of variables
# rho = 0.50
#==========================================================================================#
main.mdms.exp.cont.neg2.1 <- rsf.main(X=exp.cont.neg2.1, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cont.neg2.2 <- rsf.main(X=exp.cont.neg2.2, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cont.neg2.3 <- rsf.main(X=exp.cont.neg2.3, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cont.neg2.4 <- rsf.main(X=exp.cont.neg2.4, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cont.neg2.5 <- rsf.main(X=exp.cont.neg2.5, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cont.neg2.6 <- rsf.main(X=exp.cont.neg2.6, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cont.neg2.7 <- rsf.main(X=exp.cont.neg2.7, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cont.neg2.8 <- rsf.main(X=exp.cont.neg2.8, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cont.neg2.9 <- rsf.main(X=exp.cont.neg2.9, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cont.neg2.10 <- rsf.main(X=exp.cont.neg2.10, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                       B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)


main.mdms.exp.cat.neg2.1 <- rsf.main(X=exp.cat.neg2.1, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cat.neg2.2 <- rsf.main(X=exp.cat.neg2.2, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cat.neg2.3 <- rsf.main(X=exp.cat.neg2.3, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cat.neg2.4 <- rsf.main(X=exp.cat.neg2.4, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cat.neg2.5 <- rsf.main(X=exp.cat.neg2.5, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cat.neg2.6 <- rsf.main(X=exp.cat.neg2.6, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cat.neg2.7 <- rsf.main(X=exp.cat.neg2.7, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cat.neg2.8 <- rsf.main(X=exp.cat.neg2.8, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cat.neg2.9 <- rsf.main(X=exp.cat.neg2.9, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.exp.cat.neg2.10 <- rsf.main(X=exp.cat.neg2.10, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)


main.mdms.llv.cont.neg2.1 <- rsf.main(X=llv.cont.neg2.1, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cont.neg2.2 <- rsf.main(X=llv.cont.neg2.2, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cont.neg2.3 <- rsf.main(X=llv.cont.neg2.3, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cont.neg2.4 <- rsf.main(X=llv.cont.neg2.4, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cont.neg2.5 <- rsf.main(X=llv.cont.neg2.5, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cont.neg2.6 <- rsf.main(X=llv.cont.neg2.6, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cont.neg2.7 <- rsf.main(X=llv.cont.neg2.7, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cont.neg2.8 <- rsf.main(X=llv.cont.neg2.8, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cont.neg2.9 <- rsf.main(X=llv.cont.neg2.9, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cont.neg2.10 <- rsf.main(X=llv.cont.neg2.10, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                       B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)


main.mdms.llv.cat.neg2.1 <- rsf.main(X=llv.cat.neg2.1, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cat.neg2.2 <- rsf.main(X=llv.cat.neg2.2, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cat.neg2.3 <- rsf.main(X=llv.cat.neg2.3, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cat.neg2.4 <- rsf.main(X=llv.cat.neg2.4, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cat.neg2.5 <- rsf.main(X=llv.cat.neg2.5, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cat.neg2.6 <- rsf.main(X=llv.cat.neg2.6, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cat.neg2.7 <- rsf.main(X=llv.cat.neg2.7, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cat.neg2.8 <- rsf.main(X=llv.cat.neg2.8, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cat.neg2.9 <- rsf.main(X=llv.cat.neg2.9, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                     B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)

main.mdms.llv.cat.neg2.10 <- rsf.main(X=llv.cat.neg2.10, ntree=ntree, method="mdms", splitrule="logrank", importance="random",
                                      B=1000, ci=90, parallel=TRUE, conf=conf, verbose=TRUE, seed=12345678)


#===================================================#
# Fits a Proportional Hazards Time-To-Event Regression Model saturated with first order terms.
# Computes p-values of significance of regression coefficients of main effects in a Cox-PH model
#===================================================#
main.cph.exp.cont.neg0 <- cph.main(X=exp.cont.neg0, main.term=rownames(main.mdms.exp.cont.neg0))
main.cph.exp.cat.neg0 <- cph.main(X=exp.cat.neg0, main.term=rownames(main.mdms.exp.cat.neg0))
main.cph.llv.cont.neg0 <- cph.main(X=llv.cont.neg0, main.term=rownames(main.mdms.llv.cont.neg0))
main.cph.llv.cat.neg0 <- cph.main(X=llv.cat.neg0, main.term=rownames(main.mdms.llv.cat.neg0))

main.cph.exp.cont.neg2 <- cph.main(X=exp.cont.neg2, main.term=rownames(main.mdms.exp.cont.neg2))
main.cph.exp.cat.neg2 <- cph.main(X=exp.cat.neg2, main.term=rownames(main.mdms.exp.cat.neg2))
main.cph.llv.cont.neg2 <- cph.main(X=llv.cont.neg2, main.term=rownames(main.mdms.llv.cont.neg2))
main.cph.llv.cat.neg2 <- cph.main(X=llv.cat.neg2, main.term=rownames(main.mdms.llv.cat.neg2))

main.cph.exp.cont.neg2.1 <- cph.main(X=exp.cont.neg2.1, main.term=rownames(main.mdms.exp.cont.neg2.1))
main.cph.exp.cont.neg2.2 <- cph.main(X=exp.cont.neg2.2, main.term=rownames(main.mdms.exp.cont.neg2.2))
main.cph.exp.cont.neg2.3 <- cph.main(X=exp.cont.neg2.3, main.term=rownames(main.mdms.exp.cont.neg2.3))
main.cph.exp.cont.neg2.4 <- cph.main(X=exp.cont.neg2.4, main.term=rownames(main.mdms.exp.cont.neg2.4))
main.cph.exp.cont.neg2.5 <- cph.main(X=exp.cont.neg2.5, main.term=rownames(main.mdms.exp.cont.neg2.5))
main.cph.exp.cont.neg2.6 <- cph.main(X=exp.cont.neg2.6, main.term=rownames(main.mdms.exp.cont.neg2.6))
main.cph.exp.cont.neg2.7 <- cph.main(X=exp.cont.neg2.7, main.term=rownames(main.mdms.exp.cont.neg2.7))
main.cph.exp.cont.neg2.8 <- cph.main(X=exp.cont.neg2.8, main.term=rownames(main.mdms.exp.cont.neg2.8))
main.cph.exp.cont.neg2.9 <- cph.main(X=exp.cont.neg2.9, main.term=rownames(main.mdms.exp.cont.neg2.9))
main.cph.exp.cont.neg2.10 <- cph.main(X=exp.cont.neg2.10, main.term=rownames(main.mdms.exp.cont.neg2.10))

main.cph.exp.cat.neg2.1 <- cph.main(X=exp.cat.neg2.1, main.term=rownames(main.mdms.exp.cat.neg2.1))
main.cph.exp.cat.neg2.2 <- cph.main(X=exp.cat.neg2.2, main.term=rownames(main.mdms.exp.cat.neg2.2))
main.cph.exp.cat.neg2.3 <- cph.main(X=exp.cat.neg2.3, main.term=rownames(main.mdms.exp.cat.neg2.3))
main.cph.exp.cat.neg2.4 <- cph.main(X=exp.cat.neg2.4, main.term=rownames(main.mdms.exp.cat.neg2.4))
main.cph.exp.cat.neg2.5 <- cph.main(X=exp.cat.neg2.5, main.term=rownames(main.mdms.exp.cat.neg2.5))
main.cph.exp.cat.neg2.6 <- cph.main(X=exp.cat.neg2.6, main.term=rownames(main.mdms.exp.cat.neg2.6))
main.cph.exp.cat.neg2.7 <- cph.main(X=exp.cat.neg2.7, main.term=rownames(main.mdms.exp.cat.neg2.7))
main.cph.exp.cat.neg2.8 <- cph.main(X=exp.cat.neg2.8, main.term=rownames(main.mdms.exp.cat.neg2.8))
main.cph.exp.cat.neg2.9 <- cph.main(X=exp.cat.neg2.9, main.term=rownames(main.mdms.exp.cat.neg2.9))
main.cph.exp.cat.neg2.10 <- cph.main(X=exp.cat.neg2.10, main.term=rownames(main.mdms.exp.cat.neg2.10))

main.cph.llv.cont.neg2.1 <- cph.main(X=llv.cont.neg2.1, main.term=rownames(main.mdms.llv.cont.neg2.1))
main.cph.llv.cont.neg2.2 <- cph.main(X=llv.cont.neg2.2, main.term=rownames(main.mdms.llv.cont.neg2.2))
main.cph.llv.cont.neg2.3 <- cph.main(X=llv.cont.neg2.3, main.term=rownames(main.mdms.llv.cont.neg2.3))
main.cph.llv.cont.neg2.4 <- cph.main(X=llv.cont.neg2.4, main.term=rownames(main.mdms.llv.cont.neg2.4))
main.cph.llv.cont.neg2.5 <- cph.main(X=llv.cont.neg2.5, main.term=rownames(main.mdms.llv.cont.neg2.5))
main.cph.llv.cont.neg2.6 <- cph.main(X=llv.cont.neg2.6, main.term=rownames(main.mdms.llv.cont.neg2.6))
main.cph.llv.cont.neg2.7 <- cph.main(X=llv.cont.neg2.7, main.term=rownames(main.mdms.llv.cont.neg2.7))
main.cph.llv.cont.neg2.8 <- cph.main(X=llv.cont.neg2.8, main.term=rownames(main.mdms.llv.cont.neg2.8))
main.cph.llv.cont.neg2.9 <- cph.main(X=llv.cont.neg2.9, main.term=rownames(main.mdms.llv.cont.neg2.9))
main.cph.llv.cont.neg2.10 <- cph.main(X=llv.cont.neg2.10, main.term=rownames(main.mdms.llv.cont.neg2.10))

main.cph.llv.cat.neg2.1 <- cph.main(X=llv.cat.neg2.1, main.term=rownames(main.mdms.llv.cat.neg2.1))
main.cph.llv.cat.neg2.2 <- cph.main(X=llv.cat.neg2.2, main.term=rownames(main.mdms.llv.cat.neg2.2))
main.cph.llv.cat.neg2.3 <- cph.main(X=llv.cat.neg2.3, main.term=rownames(main.mdms.llv.cat.neg2.3))
main.cph.llv.cat.neg2.4 <- cph.main(X=llv.cat.neg2.4, main.term=rownames(main.mdms.llv.cat.neg2.4))
main.cph.llv.cat.neg2.5 <- cph.main(X=llv.cat.neg2.5, main.term=rownames(main.mdms.llv.cat.neg2.5))
main.cph.llv.cat.neg2.6 <- cph.main(X=llv.cat.neg2.6, main.term=rownames(main.mdms.llv.cat.neg2.6))
main.cph.llv.cat.neg2.7 <- cph.main(X=llv.cat.neg2.7, main.term=rownames(main.mdms.llv.cat.neg2.7))
main.cph.llv.cat.neg2.8 <- cph.main(X=llv.cat.neg2.8, main.term=rownames(main.mdms.llv.cat.neg2.8))
main.cph.llv.cat.neg2.9 <- cph.main(X=llv.cat.neg2.9, main.term=rownames(main.mdms.llv.cat.neg2.9))
main.cph.llv.cat.neg2.10 <- cph.main(X=llv.cat.neg2.10, main.term=rownames(main.mdms.llv.cat.neg2.10))

#==========================================================================================#
# RSF
# Model #5 (#pos2)
# rho = 0.50
#==========================================================================================#
formula.exp.cont.pos2 <- as.formula(paste("Surv(time=", colnames(exp.cont.pos2)[1], ", event=", colnames(exp.cont.pos2)[2], ", type=\"right\") ~ .", sep=""))
formula.exp.cat.pos2 <- as.formula(paste("Surv(time=", colnames(exp.cat.pos2)[1], ", event=", colnames(exp.cat.pos2)[2], ", type=\"right\") ~ .", sep=""))
formula.llv.cont.pos2 <- as.formula(paste("Surv(time=", colnames(llv.cont.pos2)[1], ", event=", colnames(llv.cont.pos2)[2], ", type=\"right\") ~ .", sep=""))
formula.llv.cat.pos2 <- as.formula(paste("Surv(time=", colnames(llv.cat.pos2)[1], ", event=", colnames(llv.cat.pos2)[2], ", type=\"right\") ~ .", sep=""))

rfsrc.exp.cont.pos2 <- rfsrc(formula=formula.exp.cont.pos2,
                             data=exp.cont.pos2,
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

rfsrc.exp.cat.pos2 <- rfsrc(formula=formula.exp.cat.pos2,
                             data=exp.cat.pos2,
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

rfsrc.llv.cont.pos2 <- rfsrc(formula=formula.llv.cont.pos2,
                             data=llv.cont.pos2,
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

rfsrc.llv.cat.pos2 <- rfsrc(formula=formula.llv.cat.pos2,
                             data=llv.cat.pos2,
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

#==========================================================================================
# Prediction estimates
#==========================================================================================
set.seed(12345678)
train <- sample(1:n, round(n * 4/5))
test <- setdiff(1:n, train)

#=====================
# primary call
#=====================
rfsrc.train.exp.cont.pos2 <- rfsrc(formula=formula.exp.cont.pos2,
                                   data=exp.cont.pos2[train, ],
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

rfsrc.train.exp.cat.pos2 <- rfsrc(formula=formula.exp.cat.pos2,
                                   data=exp.cat.pos2[train, ],
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

rfsrc.train.llv.cont.pos2 <- rfsrc(formula=formula.llv.cont.pos2,
                                   data=llv.cont.pos2[train, ],
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

rfsrc.train.llv.cat.pos2 <- rfsrc(formula=formula.llv.cat.pos2,
                                   data=llv.cat.pos2[train, ],
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
rfsrc.pred.exp.cont.pos2 <- predict(object=rfsrc.train.exp.cont.pos2,
                                    newdata=exp.cont.pos2[test, ],
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

rfsrc.pred.exp.cat.pos2 <- predict(object=rfsrc.train.exp.cat.pos2,
                                    newdata=exp.cat.pos2[test, ],
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

rfsrc.pred.llv.cont.pos2 <- predict(object=rfsrc.train.llv.cont.pos2,
                                    newdata=llv.cont.pos2[test, ],
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

rfsrc.pred.llv.cat.pos2 <- predict(object=rfsrc.train.llv.cat.pos2,
                                    newdata=llv.cat.pos2[test, ],
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

#==========================================================================================#
# Ranking of pairwise interactions between individual or noise variables 
# by bivariate interaction Minimal Depth of a Maximal Subtree (IMDMS)
#==========================================================================================#
int.mdms.exp.cont.neg0 <- rsf.int(X=exp.cont.neg0, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.neg1 <- rsf.int(X=exp.cont.neg1, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.neg2 <- rsf.int(X=exp.cont.neg2, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos1 <- rsf.int(X=exp.cont.pos1, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2 <- rsf.int(X=exp.cont.pos2, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos3 <- rsf.int(X=exp.cont.pos3, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos4 <- rsf.int(X=exp.cont.pos4, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos5 <- rsf.int(X=exp.cont.pos5, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos6 <- rsf.int(X=exp.cont.pos6, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos7 <- rsf.int(X=exp.cont.pos7, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)


int.mdms.exp.cat.neg0 <- rsf.int(X=exp.cat.neg0, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.neg1 <- rsf.int(X=exp.cat.neg1, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.neg2 <- rsf.int(X=exp.cat.neg2, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos1 <- rsf.int(X=exp.cat.pos1, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2 <- rsf.int(X=exp.cat.pos2, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos3 <- rsf.int(X=exp.cat.pos3, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos4 <- rsf.int(X=exp.cat.pos4, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos5 <- rsf.int(X=exp.cat.pos5, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos6 <- rsf.int(X=exp.cat.pos6, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos7 <- rsf.int(X=exp.cat.pos7, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)


int.mdms.llv.cont.neg0 <- rsf.int(X=llv.cont.neg0, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.neg1 <- rsf.int(X=llv.cont.neg1, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.neg2 <- rsf.int(X=llv.cont.neg2, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos1 <- rsf.int(X=llv.cont.pos1, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2 <- rsf.int(X=llv.cont.pos2, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos3 <- rsf.int(X=llv.cont.pos3, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos4 <- rsf.int(X=llv.cont.pos4, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos5 <- rsf.int(X=llv.cont.pos5, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos6 <- rsf.int(X=llv.cont.pos6, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos7 <- rsf.int(X=llv.cont.pos7, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)


int.mdms.llv.cat.neg0 <- rsf.int(X=llv.cat.neg0, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.neg1 <- rsf.int(X=llv.cat.neg1, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.neg2 <- rsf.int(X=llv.cat.neg2, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos1 <- rsf.int(X=llv.cat.pos1, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2 <- rsf.int(X=llv.cat.pos2, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos3 <- rsf.int(X=llv.cat.pos3, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos4 <- rsf.int(X=llv.cat.pos4, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos5 <- rsf.int(X=llv.cat.pos5, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos6 <- rsf.int(X=llv.cat.pos6, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos7 <- rsf.int(X=llv.cat.pos7, ntree=ntree, method="imdms", splitrule="logrank",
                                  importance="random", B=1000, ci=90,
                                  parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

#==========================================================================================#
# Ranking of pairwise interactions between individual or noise variables 
# by bivariate interaction Minimal Depth of a Maximal Subtree (IMDMS)
# Effect of censoring
# Regression model #2
# X1 + X2 + X1X2
#==========================================================================================#
int.mdms.exp.cont.pos2.1 <- rsf.int(X=exp.cont.pos2.1, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2.2 <- rsf.int(X=exp.cont.pos2.2, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2.3 <- rsf.int(X=exp.cont.pos2.3, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.1 <- rsf.int(X=exp.cat.pos2.1, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.2 <- rsf.int(X=exp.cat.pos2.2, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.3 <- rsf.int(X=exp.cat.pos2.3, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)


int.mdms.llv.cont.pos2.1 <- rsf.int(X=llv.cont.pos2.1, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2.2 <- rsf.int(X=llv.cont.pos2.2, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2.3 <- rsf.int(X=llv.cont.pos2.3, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.1 <- rsf.int(X=llv.cat.pos2.1, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.2 <- rsf.int(X=llv.cat.pos2.2, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.3 <- rsf.int(X=llv.cat.pos2.3, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

#==========================================================================================#
# Ranking of pairwise interactions between individual or noise variables 
# by bivariate interaction Minimal Depth of a Maximal Subtree (IMDMS)
# Effect of variable
# rho = 0.50
#==========================================================================================#
int.mdms.exp.cont.pos2.4 <- rsf.int(X=exp.cont.pos2.4, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2.5 <- rsf.int(X=exp.cont.pos2.5, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2.6 <- rsf.int(X=exp.cont.pos2.6, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2.7 <- rsf.int(X=exp.cont.pos2.7, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2.8 <- rsf.int(X=exp.cont.pos2.8, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2.9 <- rsf.int(X=exp.cont.pos2.9, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2.10 <- rsf.int(X=exp.cont.pos2.10, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2.11 <- rsf.int(X=exp.cont.pos2.11, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2.12 <- rsf.int(X=exp.cont.pos2.12, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cont.pos2.13 <- rsf.int(X=exp.cont.pos2.13, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)


int.mdms.exp.cat.pos2.4 <- rsf.int(X=exp.cat.pos2.4, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.5 <- rsf.int(X=exp.cat.pos2.5, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.6 <- rsf.int(X=exp.cat.pos2.6, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.7 <- rsf.int(X=exp.cat.pos2.7, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.8 <- rsf.int(X=exp.cat.pos2.8, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.9 <- rsf.int(X=exp.cat.pos2.9, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.10 <- rsf.int(X=exp.cat.pos2.10, ntree=ntree, method="imdms", splitrule="logrank",
                                     importance="random", B=1000, ci=90,
                                     parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.11 <- rsf.int(X=exp.cat.pos2.11, ntree=ntree, method="imdms", splitrule="logrank",
                                     importance="random", B=1000, ci=90,
                                     parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.12 <- rsf.int(X=exp.cat.pos2.12, ntree=ntree, method="imdms", splitrule="logrank",
                                     importance="random", B=1000, ci=90,
                                     parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.exp.cat.pos2.13 <- rsf.int(X=exp.cat.pos2.13, ntree=ntree, method="imdms", splitrule="logrank",
                                     importance="random", B=1000, ci=90,
                                     parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)


int.mdms.llv.cont.pos2.4 <- rsf.int(X=llv.cont.pos2.4, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2.5 <- rsf.int(X=llv.cont.pos2.5, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2.6 <- rsf.int(X=llv.cont.pos2.6, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2.7 <- rsf.int(X=llv.cont.pos2.7, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2.8 <- rsf.int(X=llv.cont.pos2.8, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2.9 <- rsf.int(X=llv.cont.pos2.9, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2.10 <- rsf.int(X=llv.cont.pos2.10, ntree=ntree, method="imdms", splitrule="logrank",
                                     importance="random", B=1000, ci=90,
                                     parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2.11 <- rsf.int(X=llv.cont.pos2.11, ntree=ntree, method="imdms", splitrule="logrank",
                                     importance="random", B=1000, ci=90,
                                     parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2.12 <- rsf.int(X=llv.cont.pos2.12, ntree=ntree, method="imdms", splitrule="logrank",
                                     importance="random", B=1000, ci=90,
                                     parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cont.pos2.13 <- rsf.int(X=llv.cont.pos2.13, ntree=ntree, method="imdms", splitrule="logrank",
                                     importance="random", B=1000, ci=90,
                                     parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)


int.mdms.llv.cat.pos2.4 <- rsf.int(X=llv.cat.pos2.4, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.5 <- rsf.int(X=llv.cat.pos2.5, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.6 <- rsf.int(X=llv.cat.pos2.6, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.7 <- rsf.int(X=llv.cat.pos2.7, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.8 <- rsf.int(X=llv.cat.pos2.8, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.9 <- rsf.int(X=llv.cat.pos2.9, ntree=ntree, method="imdms", splitrule="logrank",
                                   importance="random", B=1000, ci=90,
                                   parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.10 <- rsf.int(X=llv.cat.pos2.10, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.11 <- rsf.int(X=llv.cat.pos2.11, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.12 <- rsf.int(X=llv.cat.pos2.12, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

int.mdms.llv.cat.pos2.13 <- rsf.int(X=llv.cat.pos2.13, ntree=ntree, method="imdms", splitrule="logrank",
                                    importance="random", B=1000, ci=90,
                                    parallel=TRUE, conf=conf, verbose=FALSE, seed=seed)

#==========================================================================================#
# Fits a Proportional Hazards Time-To-Event Regression Model saturated with first and second order terms
# Computes p-values of significance of regression coefficients of pairwise interaction effects in a Cox-PH model
#==========================================================================================#
int.cph.exp.cont.neg0 <- cph.int(X=exp.cont.neg0, int.term=rownames(int.mdms.exp.cont.neg0))
int.cph.exp.cont.neg1 <- cph.int(X=exp.cont.neg1, int.term=rownames(int.mdms.exp.cont.neg1))
int.cph.exp.cont.neg2 <- cph.int(X=exp.cont.neg2, int.term=rownames(int.mdms.exp.cont.neg2))
int.cph.exp.cont.pos1 <- cph.int(X=exp.cont.pos1, int.term=rownames(int.mdms.exp.cont.pos1))
int.cph.exp.cont.pos2 <- cph.int(X=exp.cont.pos2, int.term=rownames(int.mdms.exp.cont.pos2))
int.cph.exp.cont.pos3 <- cph.int(X=exp.cont.pos3, int.term=rownames(int.mdms.exp.cont.pos3))
int.cph.exp.cont.pos4 <- cph.int(X=exp.cont.pos4, int.term=rownames(int.mdms.exp.cont.pos4))
int.cph.exp.cont.pos5 <- cph.int(X=exp.cont.pos5, int.term=rownames(int.mdms.exp.cont.pos5))
int.cph.exp.cont.pos6 <- cph.int(X=exp.cont.pos6, int.term=rownames(int.mdms.exp.cont.pos6))
int.cph.exp.cont.pos7 <- cph.int(X=exp.cont.pos7, int.term=rownames(int.mdms.exp.cont.pos7))

int.cph.exp.cont.pos2.1 <- cph.int(X=exp.cont.pos2.1, int.term=rownames(int.mdms.exp.cont.pos2.1))
int.cph.exp.cont.pos2.2 <- cph.int(X=exp.cont.pos2.2, int.term=rownames(int.mdms.exp.cont.pos2.2))
int.cph.exp.cont.pos2.3 <- cph.int(X=exp.cont.pos2.3, int.term=rownames(int.mdms.exp.cont.pos2.3))
int.cph.exp.cont.pos2.4 <- cph.int(X=exp.cont.pos2.4, int.term=rownames(int.mdms.exp.cont.pos2.4))
int.cph.exp.cont.pos2.5 <- cph.int(X=exp.cont.pos2.5, int.term=rownames(int.mdms.exp.cont.pos2.5))
int.cph.exp.cont.pos2.6 <- cph.int(X=exp.cont.pos2.6, int.term=rownames(int.mdms.exp.cont.pos2.6))
int.cph.exp.cont.pos2.7 <- cph.int(X=exp.cont.pos2.7, int.term=rownames(int.mdms.exp.cont.pos2.7))
int.cph.exp.cont.pos2.8 <- cph.int(X=exp.cont.pos2.8, int.term=rownames(int.mdms.exp.cont.pos2.8))
int.cph.exp.cont.pos2.9 <- cph.int(X=exp.cont.pos2.9, int.term=rownames(int.mdms.exp.cont.pos2.9))
int.cph.exp.cont.pos2.10 <- cph.int(X=exp.cont.pos2.10, int.term=rownames(int.mdms.exp.cont.pos2.10))
int.cph.exp.cont.pos2.11 <- cph.int(X=exp.cont.pos2.11, int.term=rownames(int.mdms.exp.cont.pos2.11))
int.cph.exp.cont.pos2.12 <- cph.int(X=exp.cont.pos2.12, int.term=rownames(int.mdms.exp.cont.pos2.12))
int.cph.exp.cont.pos2.13 <- cph.int(X=exp.cont.pos2.13, int.term=rownames(int.mdms.exp.cont.pos2.13))

int.cph.exp.cat.neg0 <- cph.int(X=exp.cat.neg0, int.term=rownames(int.mdms.exp.cat.neg0))
int.cph.exp.cat.neg1 <- cph.int(X=exp.cat.neg1, int.term=rownames(int.mdms.exp.cat.neg1))
int.cph.exp.cat.neg2 <- cph.int(X=exp.cat.neg2, int.term=rownames(int.mdms.exp.cat.neg2))
int.cph.exp.cat.pos1 <- cph.int(X=exp.cat.pos1, int.term=rownames(int.mdms.exp.cat.pos1))
int.cph.exp.cat.pos2 <- cph.int(X=exp.cat.pos2, int.term=rownames(int.mdms.exp.cat.pos2))
int.cph.exp.cat.pos3 <- cph.int(X=exp.cat.pos3, int.term=rownames(int.mdms.exp.cat.pos3))
int.cph.exp.cat.pos4 <- cph.int(X=exp.cat.pos4, int.term=rownames(int.mdms.exp.cat.pos4))
int.cph.exp.cat.pos5 <- cph.int(X=exp.cat.pos5, int.term=rownames(int.mdms.exp.cat.pos5))
int.cph.exp.cat.pos6 <- cph.int(X=exp.cat.pos6, int.term=rownames(int.mdms.exp.cat.pos6))
int.cph.exp.cat.pos7 <- cph.int(X=exp.cat.pos7, int.term=rownames(int.mdms.exp.cat.pos7))

int.cph.exp.cat.pos2.1 <- cph.int(X=exp.cat.pos2.1, int.term=rownames(int.mdms.exp.cat.pos2.1))
int.cph.exp.cat.pos2.2 <- cph.int(X=exp.cat.pos2.2, int.term=rownames(int.mdms.exp.cat.pos2.2))
int.cph.exp.cat.pos2.3 <- cph.int(X=exp.cat.pos2.3, int.term=rownames(int.mdms.exp.cat.pos2.3))
int.cph.exp.cat.pos2.4 <- cph.int(X=exp.cat.pos2.4, int.term=rownames(int.mdms.exp.cat.pos2.4))
int.cph.exp.cat.pos2.5 <- cph.int(X=exp.cat.pos2.5, int.term=rownames(int.mdms.exp.cat.pos2.5))
int.cph.exp.cat.pos2.6 <- cph.int(X=exp.cat.pos2.6, int.term=rownames(int.mdms.exp.cat.pos2.6))
int.cph.exp.cat.pos2.7 <- cph.int(X=exp.cat.pos2.7, int.term=rownames(int.mdms.exp.cat.pos2.7))
int.cph.exp.cat.pos2.8 <- cph.int(X=exp.cat.pos2.8, int.term=rownames(int.mdms.exp.cat.pos2.8))
int.cph.exp.cat.pos2.9 <- cph.int(X=exp.cat.pos2.9, int.term=rownames(int.mdms.exp.cat.pos2.9))
int.cph.exp.cat.pos2.10 <- cph.int(X=exp.cat.pos2.10, int.term=rownames(int.mdms.exp.cat.pos2.10))
int.cph.exp.cat.pos2.11 <- cph.int(X=exp.cat.pos2.11, int.term=rownames(int.mdms.exp.cat.pos2.11))
int.cph.exp.cat.pos2.12 <- cph.int(X=exp.cat.pos2.12, int.term=rownames(int.mdms.exp.cat.pos2.12))
int.cph.exp.cat.pos2.13 <- cph.int(X=exp.cat.pos2.13, int.term=rownames(int.mdms.exp.cat.pos2.13))

int.cph.llv.cont.neg0 <- cph.int(X=llv.cont.neg0, int.term=rownames(int.mdms.llv.cont.neg0))
int.cph.llv.cont.neg1 <- cph.int(X=llv.cont.neg1, int.term=rownames(int.mdms.llv.cont.neg1))
int.cph.llv.cont.neg2 <- cph.int(X=llv.cont.neg2, int.term=rownames(int.mdms.llv.cont.neg2))
int.cph.llv.cont.pos1 <- cph.int(X=llv.cont.pos1, int.term=rownames(int.mdms.llv.cont.pos1))
int.cph.llv.cont.pos2 <- cph.int(X=llv.cont.pos2, int.term=rownames(int.mdms.llv.cont.pos2))
int.cph.llv.cont.pos3 <- cph.int(X=llv.cont.pos3, int.term=rownames(int.mdms.llv.cont.pos3))
int.cph.llv.cont.pos4 <- cph.int(X=llv.cont.pos4, int.term=rownames(int.mdms.llv.cont.pos4))
int.cph.llv.cont.pos5 <- cph.int(X=llv.cont.pos5, int.term=rownames(int.mdms.llv.cont.pos5))
int.cph.llv.cont.pos6 <- cph.int(X=llv.cont.pos6, int.term=rownames(int.mdms.llv.cont.pos6))
int.cph.llv.cont.pos7 <- cph.int(X=llv.cont.pos7, int.term=rownames(int.mdms.llv.cont.pos7))

int.cph.llv.cont.pos2.1 <- cph.int(X=llv.cont.pos2.1, int.term=rownames(int.mdms.llv.cont.pos2.1))
int.cph.llv.cont.pos2.2 <- cph.int(X=llv.cont.pos2.2, int.term=rownames(int.mdms.llv.cont.pos2.2))
int.cph.llv.cont.pos2.3 <- cph.int(X=llv.cont.pos2.3, int.term=rownames(int.mdms.llv.cont.pos2.3))
int.cph.llv.cont.pos2.4 <- cph.int(X=llv.cont.pos2.4, int.term=rownames(int.mdms.llv.cont.pos2.4))
int.cph.llv.cont.pos2.5 <- cph.int(X=llv.cont.pos2.5, int.term=rownames(int.mdms.llv.cont.pos2.5))
int.cph.llv.cont.pos2.6 <- cph.int(X=llv.cont.pos2.6, int.term=rownames(int.mdms.llv.cont.pos2.6))
int.cph.llv.cont.pos2.7 <- cph.int(X=llv.cont.pos2.7, int.term=rownames(int.mdms.llv.cont.pos2.7))
int.cph.llv.cont.pos2.8 <- cph.int(X=llv.cont.pos2.8, int.term=rownames(int.mdms.llv.cont.pos2.8))
int.cph.llv.cont.pos2.9 <- cph.int(X=llv.cont.pos2.9, int.term=rownames(int.mdms.llv.cont.pos2.9))
int.cph.llv.cont.pos2.10 <- cph.int(X=llv.cont.pos2.10, int.term=rownames(int.mdms.llv.cont.pos2.10))
int.cph.llv.cont.pos2.11 <- cph.int(X=llv.cont.pos2.11, int.term=rownames(int.mdms.llv.cont.pos2.11))
int.cph.llv.cont.pos2.12 <- cph.int(X=llv.cont.pos2.12, int.term=rownames(int.mdms.llv.cont.pos2.12))
int.cph.llv.cont.pos2.13 <- cph.int(X=llv.cont.pos2.13, int.term=rownames(int.mdms.llv.cont.pos2.13))

int.cph.llv.cat.neg0 <- cph.int(X=llv.cat.neg0, int.term=rownames(int.mdms.llv.cat.neg0))
int.cph.llv.cat.neg1 <- cph.int(X=llv.cat.neg1, int.term=rownames(int.mdms.llv.cat.neg1))
int.cph.llv.cat.neg2 <- cph.int(X=llv.cat.neg2, int.term=rownames(int.mdms.llv.cat.neg2))
int.cph.llv.cat.pos1 <- cph.int(X=llv.cat.pos1, int.term=rownames(int.mdms.llv.cat.pos1))
int.cph.llv.cat.pos2 <- cph.int(X=llv.cat.pos2, int.term=rownames(int.mdms.llv.cat.pos2))
int.cph.llv.cat.pos3 <- cph.int(X=llv.cat.pos3, int.term=rownames(int.mdms.llv.cat.pos3))
int.cph.llv.cat.pos4 <- cph.int(X=llv.cat.pos4, int.term=rownames(int.mdms.llv.cat.pos4))
int.cph.llv.cat.pos5 <- cph.int(X=llv.cat.pos5, int.term=rownames(int.mdms.llv.cat.pos5))
int.cph.llv.cat.pos6 <- cph.int(X=llv.cat.pos6, int.term=rownames(int.mdms.llv.cat.pos6))
int.cph.llv.cat.pos7 <- cph.int(X=llv.cat.pos7, int.term=rownames(int.mdms.llv.cat.pos7))

int.cph.llv.cat.pos2.1 <- cph.int(X=llv.cat.pos2.1, int.term=rownames(int.mdms.llv.cat.pos2.1))
int.cph.llv.cat.pos2.2 <- cph.int(X=llv.cat.pos2.2, int.term=rownames(int.mdms.llv.cat.pos2.2))
int.cph.llv.cat.pos2.3 <- cph.int(X=llv.cat.pos2.3, int.term=rownames(int.mdms.llv.cat.pos2.3))
int.cph.llv.cat.pos2.4 <- cph.int(X=llv.cat.pos2.4, int.term=rownames(int.mdms.llv.cat.pos2.4))
int.cph.llv.cat.pos2.5 <- cph.int(X=llv.cat.pos2.5, int.term=rownames(int.mdms.llv.cat.pos2.5))
int.cph.llv.cat.pos2.6 <- cph.int(X=llv.cat.pos2.6, int.term=rownames(int.mdms.llv.cat.pos2.6))
int.cph.llv.cat.pos2.7 <- cph.int(X=llv.cat.pos2.7, int.term=rownames(int.mdms.llv.cat.pos2.7))
int.cph.llv.cat.pos2.8 <- cph.int(X=llv.cat.pos2.8, int.term=rownames(int.mdms.llv.cat.pos2.8))
int.cph.llv.cat.pos2.9 <- cph.int(X=llv.cat.pos2.9, int.term=rownames(int.mdms.llv.cat.pos2.9))
int.cph.llv.cat.pos2.10 <- cph.int(X=llv.cat.pos2.10, int.term=rownames(int.mdms.llv.cat.pos2.10))
int.cph.llv.cat.pos2.11 <- cph.int(X=llv.cat.pos2.11, int.term=rownames(int.mdms.llv.cat.pos2.11))
int.cph.llv.cat.pos2.12 <- cph.int(X=llv.cat.pos2.12, int.term=rownames(int.mdms.llv.cat.pos2.12))
int.cph.llv.cat.pos2.13 <- cph.int(X=llv.cat.pos2.13, int.term=rownames(int.mdms.llv.cat.pos2.13))

