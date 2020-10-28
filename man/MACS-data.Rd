\name{MACS}

\docType{data}

\alias{MACS}

\title{Real Dataset: MACS Cohort Study}

\description{
  Publicly available dataset from the Multicenter AIDS Cohort Study (MACS) available at (http://www.statepi.jhsph.edu/macs/macs.html).
  The dataset provides longitudinal account of viral tropism in relation to the HIV full spectrum of rates of HIV-1 disease progression (Shepherd, et al. 2008). 
  To our knowledge, this cohort provides a unique dataset with well characterized clinical information for analyzing associations between host genetic variation 
  and viral tropism as well as disease progression. Here, we determined whether copy number variation in beta-defensin and its interactions with certain polymorphisms 
  in chemokine receptors and ligand genes are associated, either alone or jointly, with clinical events in HIV-seropositive patients, 
  such as time to HIV change of tropism or time to AIDS diagnosis (See Dazard et al. (2017) for additional descriptions of the dataset and materials).
}

\usage{
    data("MACS", package="IRSF")
}

\format{
  The dataset consists of a \code{numeric} \code{data.frame} containing \eqn{n=50} complete observations (samples) 
  by rows and \eqn{p=7} covariates by columns, not including the censoring indicator and (censored) time-to-event variables. 

  The variables included in the MACS cohort study were 5 genetic variants (DEFB4/103A CNV [1-5], CCR2 SNP [190G>A], CCR5 [SNP -2459G>A, ORF], CXCL12 SNP [801G>A]) 
  and 2 non-genetic variables, taken as two additional covariates. All input variables were categorical with no more than three levels (experimental groups) each. 
  We used genetic variables with original and aggregated categories as follows: DEFB CNV [CNV = 2 or CNV > 2]; CCR2 SNP [GG or GA], CCR5 SNP [GG or GA]; 
  CCR5 ORF [WT or D32], CXCL12 SNP [GG or GA]. The first covariate was the two-level disease progression Group variable [Fast, Slow], 
  and the second was the three-level Race/Ethnicity variable [White, Hispanic, Black]. For each observation \eqn{i \in \{1,...,n\}}, we denote the \eqn{j}-th variable 
  by the \eqn{n}-dimensional vector \eqn{{\bf x}_{j} = (x_{1,j},...,x_{n,j})^{T}}, where \eqn{j \in \{1,...,p\}}. 
  Here, \eqn{p} denotes the number of variables. Hereafter, we denoted the \eqn{p=7} included variables as follows:
  \itemize{
    \item \eqn{{\bf x}_{1}}=DEFB CNV
    \item \eqn{{\bf x}_{2}}=CCR2 SNP
    \item \eqn{{\bf x}_{3}}=CCR5 SNP
    \item \eqn{{\bf x}_{4}}=CCR5 ORF
    \item \eqn{{\bf x}_{5}}=CXCL12 SNP 
    \item \eqn{{\bf x}_{6}}=Group 
    \item \eqn{{\bf x}_{7}}=Race
  }
  
  The time-to-event outcomes included in the MACS cohort study, generically denoted E, were the time-to-X4-Emergence (denoted XE) and the time-to-AIDS-Diagnosis (denoted AD), 
  whether each was observed or not during each patient's follow-up time. The corresponding event-free (EF) ("survival") probability function S(t) of time-to-event \eqn{E := XE} 
  (X4-Emergence) or \eqn{E := AD} (AIDS-Diagnosis), were called X4-Emergence-Free (\eqn{E := XEF}) or AIDS-Diagnosis-Free (\eqn{E := ADF}) probability.

  The dataset comes as a compressed Rda data file.
}

\source{
  See real data application in Dazard et al., 2017.
}

\author{
    Jean-Eudes Dazard <jean-eudes.dazard@case.edu>
    
    Maintainer: Jean-Eudes Dazard <jean-eudes.dazard@case.edu>
}

\references{
    \itemize{
        \item Dazard J-E., Ishwaran H., Mehlotra R.K., Weinberg A. and Zimmerman P.A. (2018).
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
    at (https://www.statepi.jhsph.edu/macs/macs.html) with centers at Baltimore, Chicago, Los Angeles, Pittsburgh, and the Data Coordinating Center: The Johns Hopkins University Bloomberg School of Public Health. 
    The MACS is funded primarily by the National Institute of Allergy and Infectious Diseases (NIAID), with additional co-funding from the National Cancer Institute (NCI), 
    the National Heart, Lung, and Blood Institute (NHLBI), and the National Institute on Deafness and Communication Disorders (NIDCD). MACS data collection is also supported by Johns Hopkins University CTSA. 
    This study was supported by two grants from the National Institute of Health: NIDCR P01DE019759 (Aaron Weinberg, Peter Zimmerman, Richard J. Jurevic, Mark Chance) and NCI R01CA163739 (Hemant Ishwaran). 
    The work was also partly supported by the National Science Foundation grant DMS 1148991 (Hemant Ishwaran) and the Center for AIDS Research grant P30AI036219 (Mark Chance).
}

\keyword{Real Dataset}
\keyword{Random Survival Forest}
\keyword{Interaction Detection and Modeling}
\keyword{Time-to-Event Analysis}
\keyword{Epistasis}
\keyword{Genetic Variations Interactions}

\seealso{
    \itemize{
       \item Multicenter AIDS Cohort Study (MACS) websitea at \url{https://statepi.jhsph.edu/macs/macs.html}
    }
}

\examples{
   #===================================================
   # Loading the library and its dependencies
   #===================================================
   library("IRSF")

   #===================================================
   # Help on MACS dataset
   #===================================================
   data("MACS", package="IRSF")
   ?MACS
}
