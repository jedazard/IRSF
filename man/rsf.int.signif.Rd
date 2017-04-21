\name{rsf.int.signif}
\alias{rsf.int.signif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rsf.int.signif(X, ntree, method, splitrule, importance, B, verbose, seed)
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
function (X, ntree, method, splitrule, importance, B, verbose, 
    seed) 
{
    n <- nrow(X)
    p <- ncol(X) - 2
    X.obs <- X[, -c(1, 2)]
    if (method == "mdms") {
        int.mdms.obs <- array(data = NA, dim = c(p, p, B))
        int.mdms.noise <- array(data = NA, dim = c(p, p, B))
        for (b in 1:B) {
            if (verbose) 
                cat("Iteration:", b, "\n")
            set.seed(seed = seed[b])
            X.noise <- apply(X = X.obs, MARGIN = 2, FUN = function(x) sample(x = x, 
                size = n, replace = FALSE, prob = NULL))
            colnames(X.noise) <- paste(colnames(X.obs), ".noise", 
                sep = "")
            Z.obs <- data.frame(time = X[, 1], event = X[, 2], 
                X.obs)
            Z.noise <- data.frame(time = X[, 1], event = X[, 
                2], X.noise)
            set.seed(seed = seed[b])
            bo <- sample(x = 1:n, size = n, replace = TRUE, prob = NULL)
            Z.obs.bo <- Z.obs[bo, ]
            Z.noise.bo <- Z.noise[bo, ]
            rsf.obs.bo <- rfsrc(formula = Surv(time = time, event = event, 
                type = "right") ~ ., data = Z.obs.bo, ntree = ntree, 
                bootstrap = "by.root", mtry = p, nodesize = 3, 
                splitrule = splitrule, nsplit = 0, importance = importance, 
                na.action = "na.omit", proximity = TRUE, samptype = "swr", 
                forest = TRUE, var.used = "all.trees", split.depth = "all.trees", 
                membership = TRUE, statistics = TRUE, tree.err = TRUE, 
                seed = seed[b])
            rsf.noise.bo <- rfsrc(formula = Surv(time = time, 
                event = event, type = "right") ~ ., data = Z.noise.bo, 
                ntree = ntree, bootstrap = "by.root", mtry = p, 
                nodesize = 3, splitrule = splitrule, nsplit = 0, 
                importance = importance, na.action = "na.omit", 
                proximity = TRUE, samptype = "swr", forest = TRUE, 
                var.used = "all.trees", split.depth = "all.trees", 
                membership = TRUE, statistics = TRUE, tree.err = TRUE, 
                seed = seed[b])
            int.obs.bo <- find.interaction(object = rsf.obs.bo, 
                importance = importance, method = "maxsubtree", 
                sorted = TRUE, nrep = 1, na.action = "na.omit", 
                seed = seed[b], verbose = TRUE)
            int.noise.bo <- find.interaction(object = rsf.noise.bo, 
                importance = importance, method = "maxsubtree", 
                sorted = TRUE, nrep = 1, na.action = "na.omit", 
                seed = seed[b], verbose = TRUE)
            ord.obs <- order(rownames(int.obs.bo))
            ord.noise <- order(rownames(int.noise.bo))
            int.mdms.obs[, , b] <- int.obs.bo[ord.obs, ord.obs]
            int.mdms.noise[, , b] <- int.noise.bo[ord.noise, 
                ord.noise]
        }
        dimnames(int.mdms.obs) <- list(rownames(int.obs.bo)[ord.obs], 
            colnames(int.obs.bo)[ord.obs], NULL)
        dimnames(int.mdms.noise) <- list(rownames(int.noise.bo)[ord.noise], 
            colnames(int.noise.bo)[ord.noise], NULL)
        return(list(boot.obs = int.mdms.obs, boot.noise = int.mdms.noise))
    }
    else {
        stop("Unmatched method \n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
