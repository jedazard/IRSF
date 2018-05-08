# IRSF
Interaction Random Survival Forest (IRSF): an ensemble survival tree approach to reveal variable interactions in association with time-to-events outcomes


===============
### Description
The current version is a development release. It contains R codes in folders "/R" and "/inst/doc" for the analyses and the generation of the results shown in manuscript (Dazard et al., 2017), respectively. 
Codes contain randomization, interaction modeling and prediction subroutines to be used in addition to the following R packages (http://cran.r-project.org/): [`survival`](https://CRAN.R-project.org/package=survival) for Kaplan-Meier and Cox regression modeling, [`randomForestSRC`](https://CRAN.R-project.org/package=randomForestSRC) for RSF modeling (Ishwaran and Kogalur, 2013, 2007), and [`ggRandomForests`](https://CRAN.R-project.org/package=ggRandomForests) (Ehrlinger, 2014) for Random Forrest exploration/visualization. Default parameter specifications were used for all main functions.

#### Abstract
Unraveling interactions among variables such as genetic, clinical, demographic and environmental factors is essential to understand the development of common and complex diseases. To increase the power to detect such variables interactions associated with clinical time-to-events outcomes, we borrowed established concepts from Random Survival Forest (RSF) models. We introduce a novel RSF-based pairwise interaction estimator and derive a randomization method with bootstrap confidence intervals for inferring interaction significance. Using various linear and non-linear time-to-events survival models in simulation studies, we first show the efficiency of our approach: true pairwise interaction-effects between variables are thus uncovered, while they may not be accompanied with their corresponding main-effects and often not detected by standard semi-parametric Cox regression. Moreover, using a RSF-based cross-validation scheme for generating prediction estimators, we show that informative predictors may thus be inferred. We illustrate the application of our approach in an HIV cohort study recording key host gene polymorphisms and their association with HIV change of tropism or AIDS progression. Altogether, this shows how linear or non-linear pairwise statistical interactions between variables may be uncovered in clinical studies with time-to-event outcomes of interest when the motivation is to discover important variables interactions with a predictive clinical value.

#### Key words (5)
Random Survival Forest; Interaction Detection and Modeling; Time-to-Event Analysis; Epistasis; Genetic Variations Interactions.


============
### Branches

- The default branch (master) hosts the current development release (version 1.0.2). 


===========
### License

IRSF is open source / free software, licensed under the GNU General Public License version 3 (GPLv3), 
sponsored by the [Free Software Foundation](http://www.fsf.org/). To view a copy of this license, visit 
[GNU Free Documentation License](http://www.gnu.org/licenses/gpl-3.0.html).


=============
### Downloads

(Work in progress)


================
### Requirements

IRSF (>= 1.0.2) requires R-3.0.2 (2013-09-25). It was built and tested under R version 3.5.0 (2018-04-23) and Travis CI. 

Installation has been tested on Windows, Linux, OSX and Solaris platforms. 

See Travis CI build result:
[![Build Status](https://travis-ci.org/jedazard/IRSF.png?branch=master)](https://travis-ci.org/jedazard/IRSF)


================
### Installation

* To install the most up-to-date development version (>= 1.0.2) of `IRSF` from the [GitHub](https://github.com/jedazard/IRSF) repository, 
simply run the following using devtools:

```{r}
install.packages("devtools")
library("devtools")
devtools::install_github(repo="jedazard/IRSF")
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
   + Jean-Eudes Dazard, Ph.D. [(jean-eudes.dazard@case.edu)](jean-eudes.dazard@case.edu)

Maintainers: 
   + Jean-Eudes Dazard, Ph.D. [(jean-eudes.dazard@case.edu)](jean-eudes.dazard@case.edu)

Funding/Provision/Help:   
   + This work made use of the High Performance Computing Resource in the Core Facility for Advanced Research Computing at Case Western Reserve University. 
   + We are thankful to Ms. Janet Schollenberger, Senior Project Coordinator, CAMACS, as well as Dr. Jeremy J. Martinson, Sudhir Penugonda, Shehnaz K. Hussain, Jay H. Bream, and Priya Duggal, for providing us the data related to the samples analyzed in the present study. Data in this manuscript were collected by the Multicenter AIDS Cohort Study (MACS) at (http://www.statepi.jhsph.edu/macs/macs.html) with centers at Baltimore, Chicago, Los Angeles, Pittsburgh, and the Data Coordinating Center: The Johns Hopkins University Bloomberg School of Public Health.
   + The MACS is funded primarily by the National Institute of Allergy and Infectious Diseases (NIAID), with additional co-funding from the National Cancer Institute (NCI), the National Heart, Lung, and Blood Institute (NHLBI), and the National Institute on Deafness and Communication Disorders (NIDCD). MACS data collection is also supported by Johns Hopkins University CTSA. This study was supported by two grants from the National Institute of Health: NIDCR P01DE019759 (Aaron Weinberg, Peter Zimmerman, Richard J. Jurevic, Mark Chance) and NCI R01CA163739 (Hemant Ishwaran). The work was also partly supported by the National Science Foundation grant DMS 1148991 (Hemant Ishwaran) and the Center for AIDS Research grant P30AI036219 (Mark Chance).


==============
### References

   + Dazard J-E., Ishwaran H., Mehlotra R.K., Weinberg A. and Zimmerman P.A. 
   *Ensemble Survival Tree Models to Reveal Pairwise Interactions of Variables with Time-to-Events Outcomes in Low-Dimensional Setting*. 
   [Statistical Applications in Genetics and Molecular Biology (DOI: 10.1515/sagmb-2017-0038, 2017)](https://www.degruyter.com/view/j/sagmb.2018.17.issue-1/sagmb-2017-0038/sagmb-2017-0038.xml)

   + Ishwaran, H. and Kogalur, U.B. 
   *Contributed R Package `randomForestSRC`: Random Forests for Survival, Regression and Classification (RF-SRC)*. 
   [CRAN (2013)](https://CRAN.R-project.org/package=randomForestSRC)

   + Ishwaran H. and Kogalur U.B. 
   *Random Survival Forests for R*. 
   [R News, 7(2), 25-31, (2007)](https://pdfs.semanticscholar.org/951a/84f0176076fb6786fdf43320e8b27094dcfa.pdf)

   + Ehrlinger J. 
   *Contributed R Package `ggRandomForests`: Visually Exploring Random Forests*. 
   [CRAN (2014)](https://CRAN.R-project.org/package=ggRandomForests)

   + Shepherd, J. C., et al. 
   *Emergence and Persistence of Cxcr4-Tropic Hiv-1 in a Population of Men from the Multicenter Aids Cohort Study*. 
   [J Infect Dis, 198, 1104-1112 (2008)](https://www.ncbi.nlm.nih.gov/pubmed/18783316).

