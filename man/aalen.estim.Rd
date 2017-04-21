\name{aalen.estim}
\alias{aalen.estim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
aalen.estim(obj, subset = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
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
function (obj, subset = NULL) 
{
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
            }
            else {
                r.dim <- 0
                event <- event.type <- cens <- cens <- time <- NULL
            }
            time.interest <- obj$time.interest
        }
        else {
            r.dim <- 1
            event <- event.type <- cens <- time.interest <- cens <- time <- NULL
        }
        return(list(event = event, event.type = event.type, cens = cens, 
            time.interest = time.interest, time = time, r.dim = r.dim))
    }
    event.info <- get.event.info(obj = obj, subset = subset)
    km.obj <- matrix(unlist(mclapply(X = 1:length(event.info$time.interest), 
        FUN = function(j) {
            c(sum(event.info$time >= event.info$time.interest[j], 
                na.rm = TRUE), sum(event.info$time[event.info$cens != 
                0] == event.info$time.interest[j], na.rm = TRUE))
        })), ncol = 2, byrow = TRUE)
    Y <- km.obj[, 1]
    d <- km.obj[, 2]
    r <- d/(Y + 1 * (Y == 0))
    return(exp(-cumsum(r)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
