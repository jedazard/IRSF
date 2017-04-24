\name{rsf.main.signif}

\alias{rsf.main.signif}

\title{
    Subroutine of \code{rsf.main} function
}

\description{
    Subroutine for ranking of individual and noise variables main effects by univariate Minimal Depth of a Maximal Subtree (MDMS)
}

\author{
    Jean-Eudes Dazard <jean-eudes.dazard@case.edu>
    
    Maintainer: Jean-Eudes Dazard <jean-eudes.dazard@case.edu>
}

\references{
    \itemize{
        \item Dazard J-E., Ishwaran H., Mehlotra R.K., Weinberg A. and Zimmerman P.A.
              "\emph{Ensemble Survival Tree Models to Reveal Variable Interactions in Association with Time-to-Events Outcomes}"
              Submitted (2017).
        \item Ishwaran, H. and Kogalur, U.B.
              "\emph{Random Survival Forests for R}.
              R News, 7(2), 25-31, (2007).
        \item Ishwaran, H. and Kogalur, U.B.
              "\emph{Contributed R Package \code{\link[randomForestSRC:randomForestSRC-package]{randomForestSRC}}: Random Forests for Survival, Regression and Classification (RF-SRC)}"
              CRAN (2013).
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
    \code{\link[randomForestSRC:randomForestSRC-package]{randomForestSRC}}
}
