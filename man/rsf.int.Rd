\name{rsf.int}
\alias{rsf.int}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rsf.int(X, ntree, method, splitrule = "logrank", importance = "random", B, ci = 90, parallel = FALSE, conf = NULL, verbose = TRUE, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{ntree}{
%%     ~~Describe \code{ntree} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{splitrule}{
%%     ~~Describe \code{splitrule} here~~
}
  \item{importance}{
%%     ~~Describe \code{importance} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{ci}{
%%     ~~Describe \code{ci} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{conf}{
%%     ~~Describe \code{conf} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, ntree, method, splitrule = "logrank", importance = "random", 
    B, ci = 90, parallel = FALSE, conf = NULL, verbose = TRUE, 
    seed = NULL) 
{
    p <- ncol(X) - 2
    if (!parallel) {
        if (is.null(seed)) {
            digits <- getOption("digits")
            seed <- runif(n = B, min = 1, max = 2) * 10^(digits - 
                2)
        }
        else {
            seed <- (0:(B - 1)) + seed
        }
        rsf.obj <- rsf.int.signif(X = X, ntree = ntree, method = method, 
            splitrule = splitrule, importance = importance, B = B, 
            verbose = verbose, seed = seed)
        rsf.obs.boot <- rsf.obj$boot.obs
        rsf.noise.boot <- rsf.obj$boot.noise
    }
    else {
        if (conf$type == "SOCKET") {
            cl <- parallel::makeCluster(spec = conf$spec, type = "PSOCK", 
                homogeneous = conf$homo, outfile = conf$outfile, 
                verbose = conf$verbose)
            cpus <- length(conf$spec)
        }
        else if (conf$type == "MPI") {
            cl <- parallel::makeCluster(spec = conf$spec, type = "MPI", 
                homogeneous = conf$homo, outfile = conf$outfile, 
                verbose = conf$verbose)
            cpus <- conf$spec
        }
        else {
            stop("Unrecognized cluster type\n")
        }
        parallel::clusterSetRNGStream(cl = cl, iseed = seed)
        parallel::clusterEvalQ(cl = cl, expr = library("randomForestSRC"))
        parallel::clusterEvalQ(cl = cl, expr = library("survival"))
        parallel::clusterExport(cl = cl, varlist = c("rsf.int.signif"), 
            envir = .GlobalEnv)
        rsf.cl <- parallel::clusterCall(cl = cl, fun = rsf.int.signif, 
            X = X, ntree = ntree, method = method, splitrule = splitrule, 
            importance = importance, B = ceiling(B/cpus), verbose = verbose, 
            seed = NULL)
        parallel::stopCluster(cl)
        if (method == "mdms") {
            rsf.obs.boot <- array(data = NA, dim = c(dim(rsf.cl[[1]]$boot.obs)[1], 
                dim(rsf.cl[[1]]$boot.obs)[2], 0))
            rsf.noise.boot <- array(data = NA, dim = c(dim(rsf.cl[[1]]$boot.noise)[1], 
                dim(rsf.cl[[1]]$boot.noise)[2], 0))
            for (b in 1:cpus) {
                rsf.obs.boot <- abind(rsf.obs.boot, rsf.cl[[b]]$boot.obs)
                rsf.noise.boot <- abind(rsf.noise.boot, rsf.cl[[b]]$boot.noise)
            }
        }
        else {
            stop("Unmatched method \n")
        }
    }
    theta <- (100 - ci)/200
    if (method == "mdms") {
        int.obs.mean <- apply(X = rsf.obs.boot, MARGIN = 1:2, 
            FUN = mean)
        int.noise.mean <- apply(X = rsf.noise.boot, MARGIN = 1:2, 
            FUN = mean)
        int.obs.se <- apply(X = rsf.obs.boot, MARGIN = 1:2, FUN = sd)
        int.noise.se <- apply(X = rsf.noise.boot, MARGIN = 1:2, 
            FUN = sd)
        int.obs.bpci <- apply(X = rsf.obs.boot, MARGIN = 1:2, 
            FUN = function(x) quantile(x = x, probs = c(theta, 
                1 - theta)))
        int.noise.bpci <- apply(X = rsf.noise.boot, MARGIN = 1:2, 
            FUN = function(x) quantile(x = x, probs = c(theta, 
                1 - theta)))
        int.obs.bpci <- aperm(a = int.obs.bpci, perm = c(2, 3, 
            1))
        int.noise.bpci <- aperm(a = int.noise.bpci, perm = c(2, 
            3, 1))
        vo.mean <- numeric(choose(n = p, k = 2))
        vn.mean <- numeric(choose(n = p, k = 2))
        vo.se <- numeric(choose(n = p, k = 2))
        vn.se <- numeric(choose(n = p, k = 2))
        mo.bpci <- matrix(data = NA, nrow = 4, ncol = choose(n = p, 
            k = 2))
        mn.bpci <- matrix(data = NA, nrow = 4, ncol = choose(n = p, 
            k = 2))
        k <- 1
        for (i in 2:p) {
            for (j in 1:(i - 1)) {
                vmean <- c(int.obs.mean[i, j], int.obs.mean[j, 
                  i])
                vse <- c(int.obs.se[i, j], int.obs.se[j, i])
                mbpci <- cbind(int.obs.bpci[i, j, ], int.obs.bpci[j, 
                  i, ])
                vo.mean[k] <- vmean[which.min(vmean)]
                vo.se[k] <- vse[which.min(vmean)]
                mo.bpci[, k] <- mbpci[, which.min(vmean)]
                names(vo.mean)[k] <- paste(colnames(int.obs.mean)[j], 
                  rownames(int.obs.mean)[i], sep = ":")
                vmean <- c(int.noise.mean[i, j], int.noise.mean[j, 
                  i])
                vse <- c(int.noise.se[i, j], int.noise.se[j, 
                  i])
                mbpci <- cbind(int.noise.bpci[i, j, ], int.noise.bpci[j, 
                  i, ])
                vn.mean[k] <- vmean[which.min(vmean)]
                vn.se[k] <- vse[which.min(vmean)]
                mn.bpci[, k] <- mbpci[, which.min(vmean)]
                names(vn.mean)[k] <- paste(colnames(int.noise.mean)[j], 
                  rownames(int.noise.mean)[i], sep = ":")
                k <- k + 1
            }
        }
        mat.int <- data.frame(obs.mean = vo.mean, obs.se = vo.se, 
            obs.LBCI = mo.bpci[1, ], obs.UBCI = mo.bpci[2, ], 
            noise.mean = vn.mean, noise.se = vn.se, noise.LBCI = mn.bpci[1, 
                ], noise.UBCI = mn.bpci[2, ], signif.1SE = (vo.mean + 
                vo.se < vn.mean), signif.CI = (mo.bpci[2, ] < 
                mn.bpci[1, ]))
        ord <- order(mat.int[, "obs.mean"], decreasing = F)
    }
    else {
        stop("Unmatched method \n")
    }
    mat.int <- mat.int[ord, ]
    return(mat.int)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
