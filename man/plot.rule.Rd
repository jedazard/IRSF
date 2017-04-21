\name{plot.rule}
\alias{plot.rule}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.rule(x, rule = "ci", col, lty, lwd, len = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{rule}{
%%     ~~Describe \code{rule} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{len}{
%%     ~~Describe \code{len} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, rule = "ci", col, lty, lwd, len = NULL, ...) 
{
    if (rule == "ci") {
        segments(x0 = x[, "obs.LBCI"], y0 = x[, "noise.LBCI"], 
            x1 = x[, "obs.UBCI"], y1 = x[, "noise.LBCI"], col = col, 
            lty = lty, lwd = lwd, ...)
        segments(x0 = x[, "obs.LBCI"], y0 = x[, "noise.UBCI"], 
            x1 = x[, "obs.UBCI"], y1 = x[, "noise.UBCI"], col = col, 
            lty = lty, lwd = lwd, ...)
        segments(x0 = x[, "obs.LBCI"], y0 = x[, "noise.LBCI"], 
            x1 = x[, "obs.LBCI"], y1 = x[, "noise.UBCI"], col = col, 
            lty = lty, lwd = lwd, ...)
        segments(x0 = x[, "obs.UBCI"], y0 = x[, "noise.LBCI"], 
            x1 = x[, "obs.UBCI"], y1 = x[, "noise.UBCI"], col = col, 
            lty = lty, lwd = lwd, ...)
    }
    else if (rule == "1se") {
        arrows(x0 = x[, "obs.mean"], x1 = x[, "obs.mean"] - x[, 
            "obs.se"], y0 = x[, "noise.mean"], y1 = x[, "noise.mean"], 
            length = len, angle = 90, code = 2, col = col, lty = lty, 
            lwd = lwd)
        arrows(x0 = x[, "obs.mean"], x1 = x[, "obs.mean"] + x[, 
            "obs.se"], y0 = x[, "noise.mean"], y1 = x[, "noise.mean"], 
            length = len, angle = 90, code = 2, col = col, lty = lty, 
            lwd = lwd)
    }
    else {
        warning("Not a valid rule\n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
