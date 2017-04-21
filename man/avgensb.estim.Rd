\name{avgensb.estim}
\alias{avgensb.estim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
avgensb.estim(obj, subset = NULL)
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
    if (!is.null(obj$yvar)) {
        if (is.null(subset)) {
            subset <- (1:nrow(cbind(obj$yvar)))
        }
        if (!is.null(obj$survival.oob)) {
            return(apply(obj$survival.oob[subset, , drop = FALSE], 
                2, mean, na.rm = TRUE))
        }
        else {
            return(apply(obj$survival[subset, , drop = FALSE], 
                2, mean, na.rm = TRUE))
        }
    }
    else {
        subset <- 1:obj$n
        if (!is.null(obj$survival.oob)) {
            return(apply(obj$survival.oob[subset, , drop = FALSE], 
                2, mean, na.rm = TRUE))
        }
        else {
            return(apply(obj$survival[subset, , drop = FALSE], 
                2, mean, na.rm = TRUE))
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
