\name{cph.int}

\alias{cph.int}

\title{
    Pairwise interaction effects in a Cox-PH model
}

\description{
    Fits a Proportional Hazards Time-To-Event Regression Model saturated with first and second order terms.
    Computes p-values of significance of regression coefficients of pairwise interaction effects in a Cox-PH model.
}

\usage{
    cph.int(X,
            int.term)
}

\arguments{
    \item{X}{
        \code{data.frame} or \code{numeric} \code{matrix} of input covariates.
         Dataset \code{X} assumes that:
         - all variables are in columns
         - the observed times to event and censoring variables are in the first two columns:
           \code{"stime"}: \code{numeric} \code{vector} of observed times.
           \code{"status"}: \code{numeric} \code{vector} of observed status (censoring) indicator variable.
         - each variable has a unique name, excluding the word "noise"
   }
    \item{int.term}{
    \code{vector} of \code{character} string of all possible pairs of covariates names, separated by ":".
    }
}

\value{
    \code{list} of 2 fields:
    \item{raw}{Raw p-value of covariates pairwise interaction statistics significance}
    \item{fdr}{FDR-adjusted p-value of covariates pairwise interaction statistics significance}
}

\author{
    Jean-Eudes Dazard <jean-eudes.dazard@case.edu>
    
    Maintainer: Jean-Eudes Dazard <jean-eudes.dazard@case.edu>
}

\references{
    \itemize{
        \item Dazard J-E., Ishwaran H., Mehlotra R.K., Weinberg A. and Zimmerman P.A. (2017).
              "\emph{Ensemble Survival Tree Models to Reveal Pairwise Interactions of Variables with Time-to-Events Outcomes in Low-Dimensional Setting}"
              Statistical Applications in Genetics and Molecular Biology, 17(1):20170038.
        \item Ishwaran, H. and Kogalur, U.B. (2007).
              "\emph{Random Survival Forests for R}".
              R News, 7(2):25-31.
        \item Ishwaran, H. and Kogalur, U.B. (2013).
              "\emph{Contributed R Package \pkg{randomForestSRC}: Random Forests for Survival, Regression and Classification (RF-SRC)}"
              CRAN.
    }
}

\section{Acknowledgments}{
    This work made use of the High Performance Computing Resource in the Core Facility for Advanced Research Computing at Case Western Reserve University.
    We are thankful to Ms. Janet Schollenberger, Senior Project Coordinator, CAMACS, as well as Dr. Jeremy J. Martinson, Sudhir Penugonda, Shehnaz K. Hussain, Jay H. Bream,
    and Priya Duggal, for providing us the data related to the samples analyzed in the present study. Data in this manuscript were collected by the Multicenter AIDS Cohort Study (MACS)
    at (http://www.statepi.jhsph.edu/macs/macs.html) with centers at Baltimore, Chicago, Los Angeles, Pittsburgh, and the Data Coordinating Center: The Johns Hopkins University Bloomberg School of Public Health.
    The MACS is funded primarily by the National Institute of Allergy and Infectious Diseases (NIAID), with additional co-funding from the National Cancer Institute (NCI),
    the National Heart, Lung, and Blood Institute (NHLBI), and the National Institute on Deafness and Communication Disorders (NIDCD). MACS data collection is also supported by Johns Hopkins University CTSA.
    This study was supported by two grants from the National Institute of Health: NIDCR P01DE019759 (Aaron Weinberg, Peter Zimmerman, Richard J. Jurevic, Mark Chance) and NCI R01CA163739 (Hemant Ishwaran).
    The work was also partly supported by the National Science Foundation grant DMS 1148991 (Hemant Ishwaran) and the Center for AIDS Research grant P30AI036219 (Mark Chance).
}

\keyword{Random Survival Forest}
\keyword{Interaction Detection and Modeling}
\keyword{Time-to-Event Analysis}

\seealso{
    \itemize{
       \item R package \pkg{randomForestSRC}
    }
}

\examples{
#===================================================
# Loading the library and its dependencies
#===================================================
library("IRSF")

\dontrun{
    #===================================================
    # IRSF package news
    #===================================================
    IRSF.news()

    #================================================
    # MVR package citation
    #================================================
    citation("IRSF")

    #===================================================
    # Loading of the Synthetic and Real datasets
    # Use help for descriptions
    #===================================================
    data("MACS", package="IRSF")
    ?MACS

    head(MACS)

    #===================================================
    # Synthetic dataset
    # Continuous case:
    # All variables xj, j in {1,...,p}, are iid 
    # from a multivariate uniform distribution
    # with parmeters  a=1, b=5, i.e. on [1, 5].
    # rho = 0.50
    # Regression model: X1 + X2 + X1X2
    #===================================================
    seed <- 1234567
    set.seed(seed)
    n <- 200
    p <- 5
    x <- matrix(data=runif(n=n*p, min=1, max=5),
                nrow=n, ncol=p, byrow=FALSE,
                dimnames=list(1:n, paste("X", 1:p, sep="")))

    beta <- c(rep(1,2), rep(0,p-2), 1)
    covar <- cbind(x, "X1X2"=x[,1]*x[,2])
    eta <- covar %*% beta                     # regression function

    seed <- 1234567
    set.seed(seed)
    lambda0 <- 1
    lambda <- lambda0 * exp(eta - mean(eta))  # hazards function
    tt <- rexp(n=n, rate=lambda)              # true (uncensored) event times
    tc <- runif(n=n, min=0, max=3.9)          # true (censored) event times
    stime <- pmin(tt, tc)                     # observed event times
    status <- 1 * (tt <= tc)                  # observed event indicator
    X <- data.frame(stime, status, x)

    #===================================================
    # Synthetic dataset
    # Proportional Hazards Time-To-Event Regression Model 
    # saturated with first and second order terms.
    #===================================================
    X.int.cph <- cph.int(X=X, 
                         int.term=rownames(int.mdms))
                         
    #===================================================
    # Real dataset
    #===================================================
    seed <- 1234567
    data("MACS", package="IRSF")

    X <- MACS[,c("TTX","EventX","Race","Group3",
                 "DEFB.CNV3","CCR2.SNP","CCR5.SNP2",
                 "CCR5.ORF","CXCL12.SNP2")]

    #===================================================
    # Real dataset
    # Proportional Hazards Time-To-Event Regression Model 
    # saturated with first and second order terms.
    #===================================================
    MACS.int.cph <- cph.int(X=X, 
                            int.term=rownames(int.mdms))
                         
    }
}
