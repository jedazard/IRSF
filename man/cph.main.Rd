\name{cph.main}
\alias{cph.main}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cph.main(X, main.term)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{main.term}{
%%     ~~Describe \code{main.term} here~~
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
function (X, main.term) 
{
    fmla.main <- as.formula(paste("Surv(time=", colnames(X)[1], 
        ", event=", colnames(X)[2], ", type=\"right\") ~ .", 
        sep = ""))
    P.cph.main <- numeric(p)
    names(P.cph.main) <- main.term
    for (j in 1:p) {
        Z <- X[, c(colnames(X)[1], colnames(X)[2], main.term[j])]
        coxfit <- tryCatch({
            coxph(fmla.main, data = Z, model = T, x = T, y = T)
        }, error = function(w) {
            NULL
        }, warning = function(w) {
            NULL
        })
        if (is.null(coxfit)) {
            P.cph.main[j] <- NA
        }
        else {
            P.cph.main[j] <- summary(coxfit)$coefficients[1, 
                "Pr(>|z|)", drop = T]
        }
    }
    P.cph.bh.main <- P.cph.main
    nna <- !is.na(P.cph.main)
    P.cph.bh.main[nna] <- p.adjust(p = P.cph.main[nna], method = "BH")
    P.cph.bh.main[!nna] <- NA
    return(list(raw = P.cph.main, fdr = P.cph.bh.main))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
