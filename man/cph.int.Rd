\name{cph.int}
\alias{cph.int}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cph.int(X, int.term)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{int.term}{
%%     ~~Describe \code{int.term} here~~
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
function (X, int.term) 
{
    p <- dim(X)[2] - 2
    int.term.list <- strsplit(x = int.term, split = ":")
    fmla.int <- as.formula(paste("Surv(time=", colnames(X)[1], 
        ", event=", colnames(X)[2], ", type=\"right\") ~ . + (.)^2", 
        sep = ""))
    P.cph.int <- numeric(p * (p - 1)/2)
    names(P.cph.int) <- int.term
    for (j in 1:(p * (p - 1)/2)) {
        Z <- X[, c(colnames(X)[1], colnames(X)[2], int.term.list[[j]])]
        coxfit <- tryCatch({
            coxph(fmla.int, data = Z, model = T, x = T, y = T)
        }, error = function(w) {
            NULL
        }, warning = function(w) {
            NULL
        })
        if (is.null(coxfit)) {
            P.cph.int[j] <- NA
        }
        else {
            P.cph.int[j] <- summary(coxfit)$coefficients[3, "Pr(>|z|)", 
                drop = T]
        }
    }
    P.cph.bh.int <- P.cph.int
    nna <- !is.na(P.cph.int)
    P.cph.bh.int[nna] <- p.adjust(p = P.cph.int[nna], method = "BH")
    P.cph.bh.int[!nna] <- NA
    return(list(raw = P.cph.int, fdr = P.cph.bh.int))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
