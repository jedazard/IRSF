\name{rsf.rank}
\alias{rsf.rank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rsf.rank(X, ntree, method, splitrule = "logrank", importance = "random", B, ci = 90, parallel = FALSE, conf = NULL, verbose = TRUE, seed = NULL)
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
    if (!parallel) {
        if (is.null(seed)) {
            digits <- getOption("digits")
            seed <- runif(n = B, min = 1, max = 2) * 10^(digits - 
                2)
        }
        else {
            seed <- (0:(B - 1)) + seed
        }
        rsf.obj <- rsf.rank.signif(X = X, ntree = ntree, method = method, 
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
        parallel::clusterExport(cl = cl, varlist = c("rsf.rank.signif"), 
            envir = .GlobalEnv)
        rsf.cl <- parallel::clusterCall(cl = cl, fun = rsf.rank.signif, 
            X = X, ntree = ntree, method = method, splitrule = splitrule, 
            importance = importance, B = ceiling(B/cpus), verbose = verbose, 
            seed = NULL)
        parallel::stopCluster(cl)
        rsf.obs.boot <- matrix(data = NA, nrow = 0, ncol = ncol(rsf.cl[[1]]$boot.obs))
        rsf.noise.boot <- matrix(data = NA, nrow = 0, ncol = ncol(rsf.cl[[1]]$boot.noise))
        for (b in 1:cpus) {
            rsf.obs.boot <- rbind(rsf.obs.boot, rsf.cl[[b]]$boot.obs)
            rsf.noise.boot <- rbind(rsf.noise.boot, rsf.cl[[b]]$boot.noise)
        }
    }
    theta <- (100 - ci)/200
    ranks.obs.mean <- apply(X = rsf.obs.boot, MARGIN = 2, FUN = mean)
    ranks.noise.mean <- apply(X = rsf.noise.boot, MARGIN = 2, 
        FUN = mean)
    ranks.obs.se <- apply(X = rsf.obs.boot, MARGIN = 2, FUN = sd)
    ranks.noise.se <- apply(X = rsf.noise.boot, MARGIN = 2, FUN = sd)
    ranks.obs.bpci <- apply(X = rsf.obs.boot, MARGIN = 2, FUN = function(x) quantile(x = x, 
        probs = c(theta, 1 - theta)))
    ranks.noise.bpci <- apply(X = rsf.noise.boot, MARGIN = 2, 
        FUN = function(x) quantile(x = x, probs = c(theta, 1 - 
            theta)))
    names(ranks.noise.mean) <- gsub(pattern = ".noise", replacement = "", 
        x = names(ranks.noise.mean), ignore.case = F, fixed = F)
    w <- pmatch(x = names(ranks.obs.mean), table = names(ranks.noise.mean))
    ranks.noise.mean <- ranks.noise.mean[w]
    ranks.noise.se <- ranks.noise.se[w]
    ranks.noise.bpci <- ranks.noise.bpci[, w]
    if (method == "mdms") {
        mat.ranks <- data.frame(obs.mean = ranks.obs.mean, obs.se = ranks.obs.se, 
            obs.LBCI = ranks.obs.bpci[1, ], obs.UBCI = ranks.obs.bpci[2, 
                ], noise.mean = ranks.noise.mean, noise.se = ranks.noise.se, 
            noise.LBCI = ranks.noise.bpci[1, ], noise.UBCI = ranks.noise.bpci[2, 
                ], signif.1SE = (ranks.obs.mean + ranks.obs.se < 
                ranks.noise.mean), signif.CI = (ranks.obs.bpci[2, 
                ] < ranks.noise.bpci[1, ]))
        ord <- order(mat.ranks[, "obs.mean"], decreasing = F)
    }
    else {
        stop("Unmatched method \n")
    }
    mat.ranks <- mat.ranks[ord, ]
    return(mat.ranks)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
