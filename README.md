# IRSF
Interaction Random Survival Forest (IRSF): an ensemble survival tree approach to reveal variable interactions in association with time-to-events outcomes


===============
### Description
Builds ensemble survival tree models to reveal variable interactions when the response is a time-to-events outcome. 
Codes contain randomization, interaction modeling, and prediction subroutines to be used in addition to the following 
R packages: [`survival`](https://CRAN.R-project.org/package=survival) for Kaplan-Meier and Cox regression modeling, 
[`randomForestSRC`](https://CRAN.R-project.org/package=randomForestSRC) (Ishwaran and Kogalur, 2013, 2007) for RSF modeling, 
and optionally [`ggRandomForests`](https://CRAN.R-project.org/package=ggRandomForests) (Ehrlinger, 2014) for Random Forest 
exploration/visualization. The current version contains additional R codes in folder "/inst/doc" for the analysis and generation 
of results shown in the corresponding article (Dazard et al., 2018).

#### Abstract
Unraveling interactions among variables such as genetic, clinical, demographic and environmental factors is essential to understand the development of common and complex diseases. To increase the power to detect such variables interactions associated with clinical time-to-events outcomes, we borrowed established concepts from Random Survival Forest (RSF) models. We introduce a novel RSF-based pairwise interaction estimator and derive a randomization method with bootstrap confidence intervals for inferring interaction significance. Using various linear and non-linear time-to-events survival models in simulation studies, we first show the efficiency of our approach: true pairwise interaction-effects between variables are thus uncovered, while they may not be accompanied with their corresponding main-effects and often not detected by standard semi-parametric Cox regression. Moreover, using a RSF-based cross-validation scheme for generating prediction estimators, we show that informative predictors may thus be inferred. We illustrate the application of our approach in an HIV cohort study recording key host gene polymorphisms and their association with HIV change of tropism or AIDS progression. Altogether, this shows how linear or non-linear pairwise statistical interactions between variables may be uncovered in clinical studies with time-to-event outcomes of interest when the motivation is to discover important variables interactions with a predictive clinical value.

#### Key words (5)
Random Survival Forest; Interaction Detection and Modeling; Time-to-Event Analysis; Epistasis; Genetic Variations Interactions.


============
### Branches

This branch (master) is the  default one, that hosts the current development release (version 1.0.3).


===========
### License

IRSF is open source / free software, licensed under the GNU General Public License version 3 (GPLv3), 
sponsored by the [Free Software Foundation](https://www.fsf.org/). To view a copy of this license, visit 
[GNU Free Documentation License](https://www.gnu.org/licenses/gpl-3.0.html).


=============
### Downloads

CRAN downloads since October 1, 2012, 
the month the [RStudio CRAN mirror](http://cran-logs.rstudio.com/) 
started publishing logs:
[![](https://cranlogs.r-pkg.org/badges/grand-total/IRSF)](https://CRAN.R-project.org/package=IRSF)

CRAN downloads in the last month:
[![](https://cranlogs.r-pkg.org/badges/last-month/IRSF)](https://CRAN.R-project.org/package=IRSF)

CRAN downloads in the last week:
[![](https://cranlogs.r-pkg.org/badges/last-week/IRSF)](https://CRAN.R-project.org/package=IRSF)


================
### Requirements

IRSF (>= 1.0.3) requires R-3.5.0 (2018-04-23). It was built and tested under R version 4.0.3 (2020-10-10) and Travis CI. 

Installation has been tested on Windows, Linux, OSX and Solaris platforms. 

See Travis CI build result:
[![Build Status](https://app.travis-ci.com/jedazard/IRSF.svg)](https://app.travis-ci.com/jedazard/IRSF)

See CRAN checks:
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/IRSF)](https://cran.r-project.org/web/checks/check_results_IRSF.html)


================
### Installation

* To install the stable version of `IRSF`, simply download and install the current version (1.0.3) from the [CRAN](https://CRAN.R-project.org/package=IRSF) 
repository:

```{r}
install.packages("IRSF")
```

* Alternatively, you can install the most up-to-date development version (>= 1.0.3) of `IRSF` from the [GitHub](https://github.com/jedazard/IRSF) repository:

```{r}
install.packages("devtools")
library("devtools")
devtools::install_github("jedazard/IRSF")
```

=========
### Usage

* To load the IRSF library in an R session and start using it:

```{r}
library("IRSF")
```

* Check on how to cite the package with the R command:

```{r}
citation("IRSF")
```

etc...


===================
### Acknowledgments

Authors: 
   + Jean-Eudes Dazard, Ph.D. <jean-eudes.dazard@case.edu>

Maintainers: 
   + Jean-Eudes Dazard, Ph.D. <jean-eudes.dazard@case.edu>

Funding/Provision/Help:
   + This work made use of the High Performance Computing Resource in the Core Facility for Advanced Research Computing at Case Western Reserve University. 
   + We are thankful to Ms. Janet Schollenberger, Senior Project Coordinator, CAMACS, as well as Dr. Jeremy J. Martinson, Sudhir Penugonda, Shehnaz K. Hussain, Jay H. Bream, and Priya Duggal, for providing us the data related to the samples analyzed in the present study. Data in this manuscript were collected by the Multicenter AIDS Cohort Study (MACS) at (https://statepi.jhsph.edu/macs/macs.html) with centers at Baltimore, Chicago, Los Angeles, Pittsburgh, and the Data Coordinating Center: The Johns Hopkins University Bloomberg School of Public Health.
   + The MACS is funded primarily by the National Institute of Allergy and Infectious Diseases (NIAID), with additional co-funding from the National Cancer Institute (NCI), the National Heart, Lung, and Blood Institute (NHLBI), and the National Institute on Deafness and Communication Disorders (NIDCD). MACS data collection is also supported by Johns Hopkins University CTSA. This study was supported by two grants from the National Institute of Health: NIDCR P01DE019759 (Aaron Weinberg, Peter Zimmerman, Richard J. Jurevic, Mark Chance) and NCI R01CA163739 (Hemant Ishwaran). The work was also partly supported by the National Science Foundation grant DMS 1148991 (Hemant Ishwaran) and the Center for AIDS Research grant P30AI036219 (Mark Chance).


==============
### References

   + Dazard J-E., Ishwaran H., Mehlotra R.K., Weinberg A. and Zimmerman P.A. 
   *Ensemble Survival Tree Models to Reveal Pairwise Interactions of Variables with Time-to-Events Outcomes in Low-Dimensional Setting*. 
   [Statistical Applications in Genetics and Molecular Biology (2018), 17(1):20170038](https://doi.org/10.1515/sagmb-2017-0038).

   + Ishwaran, H. and Kogalur, U.B. 
   *Contributed R Package `randomForestSRC`: Random Forests for Survival, Regression and Classification (RF-SRC)*. 
   [CRAN (2013)](https://CRAN.R-project.org/package=randomForestSRC)

   + Ishwaran H. and Kogalur U.B. 
   *Random Survival Forests for R*. 
   [R News, 7(2), 25-31, (2007)](https://www.semanticscholar.org/paper/Random-survival-forests-Ishwaran-Kogalur/9ee2d6a8de063e2621eebc620b9d9d3d8a380374?p2df)

   + Ehrlinger J. 
   *Contributed R Package `ggRandomForests`: Visually Exploring Random Forests*. 
   [CRAN (2014)](https://CRAN.R-project.org/package=ggRandomForests)

   + Shepherd, J. C., et al. 
   *Emergence and Persistence of Cxcr4-Tropic Hiv-1 in a Population of Men from the Multicenter Aids Cohort Study*. 
   [J Infect Dis, 198, 1104-1112 (2008)](https://pubmed.ncbi.nlm.nih.gov/18783316/).

