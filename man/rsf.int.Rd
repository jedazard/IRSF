\name{rsf.int}

\alias{rsf.int}

\title{
    Bivariate Interaction Minimal Depth of a Maximal Subtree (IMDMS)
}

\description{
    Ranking of pairwise interactions between individual or noise variables by bivariate interaction Minimal Depth of a Maximal Subtree (IMDMS)
}

\usage{
    rsf.int(X,
            ntree = 1000,
            method = "imdms",
            splitrule = "logrank",
            importance = "random",
            B, ci = 90,
            parallel = FALSE,
            conf = NULL,
            verbose = TRUE,
            seed = NULL)
}

\arguments{
    \item{X}{
        \code{data.frame} or \code{numeric} \code{matrix} of input covariates.
         Dataset \code{X} assumes that:
         - all variables are in columns
         - the observed times to event and censoring variables are in the first two columns:
           "stime": \code{numeric} \code{vector} of observed times.
           "status": \code{numeric} \code{vector} of observed status (censoring) indicator variable.
         - each variable has a unique name, excluding the word "noise"
    }
    \item{ntree}{
        Number of trees in the forest. Defaults to 1000.
    }
    \item{method}{
        Method for ranking of interactions between pairs of individual and noise variables.
        \code{character} string \code{"imdms"} (default) that stands for interaction minimal depth of a maximal subtree (IMDMS).
    }
    \item{splitrule}{
        Splitting rule used to grow trees. For time-to-event analysis, use \code{"logrank"} (default),
        which implements log-rank splitting (Segal, 1988; Leblanc and Crowley, 1993).
    }
    \item{importance}{
        Method for computing variable importance. Defaults to Character string \code{"random"}. See details below.
    }
    \item{B}{
        Postitive \code{integer} of the number of replications of the cross-validation procedure.
    }
    \item{ci}{
        Confidence Interval for inferences of individual and noise variables.
        \code{numeric} scalar between 50 and 100. Defaults to 90.
    }
    \item{parallel}{
        \code{logical}. Is parallel computing to be performed? Defaults to \code{FALSE}.
    }
    \item{conf}{
        \code{list} of 5 fields containing the parameters values needed for creating the parallel backend (cluster configuration).
        See details below for usage. Optional, defaults to \code{NULL}, but all fields are required if used:
        \itemize{
            \item \code{type} : \code{character} \code{vector} specifying the cluster type ("SOCKET", "MPI").
            \item \code{spec} : A specification (\code{character} \code{vector} or \code{integer} scalar) appropriate to the type of cluster.
            \item \code{homogeneous} : \code{logical} scalar to be set to \code{FALSE} for inhomogeneous clusters.
            \item \code{verbose} : \code{logical} scalar to be set to \code{FALSE} for quiet mode.
            \item \code{outfile} : \code{character} \code{vector} of an output log file name to direct the stdout and stderr
                  connection output from the workernodes. "" indicates no redirection.
        }
    }
    \item{verbose}{
        \code{logical} scalar. Is the output to be verbose? Optional, defaults to \code{TRUE}.
    }
    \item{seed}{
        Positive \code{integer} scalar of the user seed to reproduce the results. Defaults to \code{NULL}.
    }
}

\details{
    The option \code{importance} allows several ways to calculate Variable Importance (VIMP).
    The default \code{"permute"} returns Breiman-Cutler permutation VIMP as described in Breiman (2001).
    For each tree, the prediction error on the out-of-bag (OOB) data is recorded.
    Then for a given variable \code{x}, OOB cases are randomly permuted in \code{x} and the prediction error is recorded.
    The VIMP for \code{x} is defined as the difference between the perturbed and unperturbed error rate,
    averaged over all trees. If \code{"random"} is used, then \code{x} is not permuted,
    but rather an OOB case is assigned a daughter node randomly whenever a split on \code{x} is encountered in the in-bag tree.
    If \code{"anti"} is used, then \code{x} is assigned to the opposite node whenever a split on \code{x} is encountered in the in-bag tree.

    The function \code{rsf.int} relies on the R package \pkg{parallel} to create a parallel backend within an R session, enabling access to a cluster
    of compute cores and/or nodes on a local and/or remote machine(s) and scaling-up with the number of CPU cores available and efficient parallel
    execution. To run a procedure in parallel (with parallel RNG), argument \code{parallel} is to be set to \code{TRUE} and argument \code{conf}
    is to be specified (i.e. non \code{NULL}). Argument \code{conf} uses the options described in function \code{makeCluster} of the R packages
    \pkg{parallel} and \pkg{snow}. \pkg{IRSF} supports two types of communication mechanisms between master and worker processes:
    'Socket' or  'Message-Passing Interface' ('MPI'). In \pkg{IRSF}, parallel 'Socket' clusters use sockets communication mechanisms only
    (no forking) and are therefore available on all platforms, including Windows, while parallel 'MPI' clusters use high-speed interconnects
    mechanism in networks of computers (with distributed memory) and are therefore available only in these architectures. A parallel 'MPI'
    cluster also requires R package \pkg{Rmpi} to be installed. Value \code{type} is used to setup a cluster of type 'Socket' ("SOCKET")
    or 'MPI' ("MPI"), respectively. Depending on this type, values of \code{spec} are to be used alternatively:
    \itemize{
        \item For 'Socket' clusters (\code{conf$type="SOCKET"}), \code{spec} should be a \code{character} \code{vector} naming the hosts on which
              to run the job; it can default to a unique local machine, in which case, one may use the unique host name "localhost".
              Each host name can potentially be repeated to the number of CPU cores available on the local machine.
              It can also be an \code{integer} scalar specifying the number of processes to spawn on the local machine;
              or a list of machine specifications if you have ssh installed (a character value named host specifying the name or address of the host to use).
        \item For 'MPI' clusters (\code{conf$type="MPI"}), \code{spec} should be an \code{integer} scalar
              specifying the total number of processes to be spawned across the network of available nodes, counting the workernodes and masternode.
    }
    The actual creation of the cluster, its initialization, and closing are all done internally. For more details,
    see the reference manual of R package \pkg{snow} and examples below.

    When random number generation is needed, the creation of separate streams of parallel RNG per node is done internally by
    distributing the stream states to the nodes. For more details, see the vignette of R package \pkg{parallel}.
    The use of a seed allows to reproduce the results within the same type of session: the same seed will reproduce the same results within
    a non-parallel session or within a parallel session, but it will not necessarily give the exact same results (up to sampling variability)
    between a non-parallelized and parallelized session due to the difference of management of the seed between the two (see parallel RNG and
    value of returned seed below).
}

\value{
    \code{data.frame} containing the following columns:
    \itemize{ 
        \item "obs.mean"    observed mean of covariates pairwise interaction statistics
        \item "obs.se"      observed standard error of covariates pairwise interaction statistics
        \item "obs.LBCI"    observed Lower Bound Confidence Interval of covariates pairwise interaction ranstatisticsks
        \item "obs.UBCI"    observed Upper Bound Confidence Interval of covariates pairwise interaction statistics
        \item "noise.mean"  observed mean of noise covariates statistics
        \item "noise.se"    observed standard error of noise covariates pairwise interaction statistics
        \item "noise.LBCI"  observed Lower Bound Confidence Interval of noise covariates pairwise interaction statistics
        \item "noise.UBCI"  observed Upper Bound Confidence Interval of noise covariates pairwise interaction statistics
        \item "signif.1SE"  calls of covariates pairwise interaction statistics significance using the 1SE rule
        \item "signif.CI"   calls of covariates pairwise interaction statistics significance using the CI rule at \code{ci}\% confidence level
    }
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

\examples{
   #===================================================
   # Loading the library and its dependencies
   #===================================================
   library("IRSF")

   #==========================================================================================#
   # Continuous case:
   # All variables xj, j in {1,...,p}, are iid from a multivariate uniform distribution
   # with parmeters  a=1, b=5, i.e. on [1, 5].
   # rho = 0.50
   # Regression model: X1 + X2 + X1X2
   #==========================================================================================#
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

   int.mdms <- rsf.int(X=X, 
                       ntree=1000, 
                       method="imdms", 
                       splitrule="logrank",
                       importance="random", 
                       B=5, 
                       ci=90,
                       parallel=FALSE, 
                       conf=NULL, 
                       verbose=FALSE, 
                       seed=seed)

   \dontrun{
      #===================================================
      # Examples of parallel backend parametrization
      #===================================================
      # Example #1 - Quad core PC
      # Running WINDOWS with SOCKET communication
      #===================================================
      cpus <- detectCores(logical = TRUE)
      conf <- list("spec" = rep("localhost", cpus),
                   "type" = "SOCK",
                   "homo" = TRUE,
                   "verbose" = TRUE,
                   "outfile" = "")
      #===================================================
      # Example #2 - Master node + 3 Worker nodes cluster
      # Running LINUX with SOCKET communication
      # All nodes equipped with identical setups of
      # multicores (8 core CPUs per machine for a total of 32)
      #===================================================
      masterhost <- Sys.getenv("HOSTNAME")
      slavehosts <- c("compute-0-0", "compute-0-1", "compute-0-2")
      nodes <- length(slavehosts) + 1
      cpus <- 8
      conf <- list("spec" = c(rep(masterhost, cpus),
                              rep(slavehosts, cpus)),
                   "type" = "SOCK",
                   "homo" = TRUE,
                   "verbose" = TRUE,
                   "outfile" = "")
      #===================================================
      # Example #3 - Multinode of multicore per node cluster
      # Running LINUX with SLURM scheduler and MPI communication
      # Below, variable 'cpus' is the total number
      # of requested core CPUs, which is specified from
      # within a SLURM script.
      #===================================================
      cpus <- as.numeric(Sys.getenv("SLURM_NTASKS"))
      if (require("Rmpi")) {
         print("Rmpi is loaded correctly \n")
      } else {
         stop("Rmpi must be installed first \n")
      }
      conf <- list("spec" = cpus,
                   "type" = "MPI",
                   "homo" = TRUE,
                   "verbose" = TRUE,
                   "outfile" = "")
   }
}
