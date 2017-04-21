\name{IRSF-package}

\alias{IRSF}

\docType{package}

\title{
\packageTitle{IRSF}
}

\description{
\packageDescription{IRSF}
}

\details{
Unraveling interactions among variables such as genetic, clinical, demographic and environmental factors 
is essential to understand the development of common and complex diseases. To increase the power to detect 
such variables interactions associated with clinical time-to-events outcomes, we borrowed established 
concepts from Random Survival Forest (RSF) models. We introduce a novel RSF-based pairwise interaction 
estimator and derive a randomization method with bootstrap confidence intervals for inferring interaction significance. 
Using various linear and non-linear time-to-events survival models in simulation studies, we first show the 
efficiency of our approach: true pairwise interaction-effects between variables are thus uncovered, 
while they may not be accompanied with their corresponding main-effects and often not detected 
by standard semi-parametric Cox regression. Moreover, using a RSF-based cross-validation scheme 
for generating prediction estimators, we show that informative predictors may thus be inferred. 
We illustrate the application of our approach in an HIV cohort study recording key host gene polymorphisms 
and their association with HIV change of tropism or AIDS progression. Altogether, this shows how linear 
or non-linear pairwise statistical interactions between variables may be uncovered in clinical studies 
with time-to-event outcomes of interest when the motivation is to discover important variables interactions
with a predictive clinical value.

The DESCRIPTION file:
\packageDESCRIPTION{IRSF}
\packageIndices{IRSF}
~~ An overview of how to use the package, including the most important functions ~~
}

\author{
\packageAuthor{IRSF}
Maintainer: \packageMaintainer{IRSF}
}

\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
~~ documentation directory ~~

\keyword{Random Survival Forest} 
\keyword{Interaction Detection and Modeling} 
\keyword{Time-to-Event Analysis} 
\keyword{Epistasis} 
\keyword{Genetic Variations Interactions}

\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}

