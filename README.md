# IRSF
Interaction Random Survival Forest (IRSF): Ensemble Survival Tree Models to Reveal Variable Interactions in Association with Time-to-Events Outcomes


==============
### License

PRIMsrc is open source / free software, licensed under the GNU General Public License version 3 (GPLv3), 
sponsored by the [Free Software Foundation](http://www.fsf.org/). To view a copy of this license, visit 
[GNU Free Documentation License](http://www.gnu.org/licenses/gpl-3.0.html).


==============
### Description
The current version is a development release. It contains R codes for the analyses and the generation of the results shown in manuscript (Dazard et al., 2017). Codes contain randomization, interaction modeling and prediction subroutines to be used in addition to the following R packages (http://cran.r-project.org/): [`survival`](https://CRAN.R-project.org/package=survival) for Kaplan-Meier and Cox regression modeling, [`NADA`](https://CRAN.R-project.org/package=NADA) for correlation analysis in the presence of censoring, [`randomForestSRC`](https://CRAN.R-project.org/package=randomForestSRC) for RSF modeling (Ishwaran and Kogalur, 2013, 2007), and [`ggRandomForests`](https://CRAN.R-project.org/package=ggRandomForests) (Ehrlinger, 2014) for Random Forrest exploration/visualization. Default parameter specifications were used for all main functions.

==============
### Abstract
Unraveling interactions among variables such as genetic, clinical, demographic and environmental factors is essential to understand the development of common and complex diseases. To increase the power to detect such variables interactions associated with clinical time-to-events outcomes, we borrowed established concepts from Random Survival Forest (RSF) models. We introduce a novel RSF-based pairwise interaction estimator and derive a randomization method with bootstrap confidence intervals for inferring interaction significance. Using various linear and non-linear time-to-events survival models in simulation studies, we first show the efficiency of our approach: true pairwise interaction-effects between variables are thus uncovered, while they may not be accompanied with their corresponding main-effects and often not detected by standard semi-parametric Cox regression. Moreover, using a RSF-based cross-validation scheme for generating prediction estimators, we show that informative predictors may thus be inferred. We illustrate the application of our approach in an HIV cohort study recording key host gene polymorphisms and their association with HIV change of tropism or AIDS progression. Altogether, this shows how linear or non-linear pairwise statistical interactions between variables may be uncovered in clinical studies with time-to-event outcomes of interest when the motivation is to discover important variables interactions with a predictive clinical value.

##### Key words (5)
Random Survival Forest; Interaction Detection and Modeling; Time-to-Event Analysis; Epistasis; Genetic Variations Interactions.


==============
### Real Dataset
##### Background
For our objectives, we utilized samples from the previously published MACS cohort study (http://www.statepi.jhsph.edu/macs/macs.html), which provides longitudinal account of viral tropism in relation to the HIV full spectrum of rates of HIV-1 disease progression (Shepherd, et al. 2008). 
To our knowledge, this cohort provides a unique dataset with well characterized clinical information for analyzing associations between host genetic variation and viral tropism as well as disease progression. 
Here, we determined whether copy number variation in beta-defensin and its interactions with certain polymorphisms in chemokine receptors and ligand genes are associated, either alone or jointly, with clinical events in HIV-seropositive patients, such as time to HIV change of tropism or time to AIDS diagnosis. 
Additional descriptions of the dataset and materials used are provided in the Supplemental Materials.

##### Variables and Outcomes 
The variables included in the MACS cohort study were 5 genetic variants (DEFB4/103A CNV [1-5], CCR2 SNP [190G>A], CCR5 [SNP -2459G>A, ORF], CXCL12 SNP [801G>A]) and 2 non-genetic variables, taken as two additional covariates. 
All input variables were categorical with no more than three levels (experimental groups) each. We used genetic variables with original and aggregated categories as follows: DEFB CNV [CNV = 2 or CNV > 2]; CCR2 SNP [GG or GA], CCR5 SNP [GG or GA]; CCR5 ORF [WT or D32], CXCL12 SNP [GG or GA]. 
The first covariate was the two-level disease progression Group variable ["Fast", "Slow"], and the second was the three-level Race/Ethnicity variable [White, Hispanic, Black]. 
For each observation _i_ \in {1,...,_n_}, we denote the _j_-th variable by the _n_-dimensional vector **x**<sub>_j_</sub> = (x<sub>1,_j_</sub>,...,x<sub>_n_,_j_</sub>)<sup>_T_</sup>, where _j_ \in {1,...,_p_}. 
Here, _p_ denotes the number of variables. Hereafter, we denoted the _p_= 7 included variables as follows: **x**<sub>1</sub>=DEFB CNV, **x**<sub>2</sub>=CCR2 SNP, **x**<sub>3</sub>=CCR5 SNP, **x**<sub>4</sub>=CCR5 ORF, **x**<sub>5</sub>=CXCL12 SNP, **x**<sub>6</sub>=Group, **x**<sub>7</sub>=Race.

The time-to-event outcomes included in the MACS cohort study, generically denoted E, were the time-to-X4-Emergence (denoted XE) and the time-to-AIDS-Diagnosis (denoted AD), whether each was observed or not during each patient’s follow-up time. The corresponding event-free (EF) (“survival”) probability function S(t) of time-to-event E := XE (X4-Emergence) or E := AD (AIDS-Diagnosis),were called X4-Emergence-Free (E := XEF) or AIDS-Diagnosis-Free (E := ADF) probability.

==============
### References
- Dazard J-E., Ishwaran H., Mehlotra R.K., Weinberg A. and Zimmerman P.A. Ensemble Survival Tree Models to Reveal Variable Interactions in Association with Time-to-Events Outcomes. 
[submitted (2017)]

- Ishwaran, H. and Kogalur, U.B. Contributed R Package `randomForestSRC`: Random Forests for Survival, Regression and Classification (RF-SRC).
[CRAN (2013)](https://CRAN.R-project.org/package=randomForestSRC)

- Ishwaran H. and Kogalur U.B. Random Survival Forests for R. 
[R News, 7(2), 25-31, (2007)](https://pdfs.semanticscholar.org/951a/84f0176076fb6786fdf43320e8b27094dcfa.pdf)

- Ehrlinger J. Contributed R Package `ggRandomForests`: Visually Exploring Random Forests
[CRAN (2014)](https://CRAN.R-project.org/package=ggRandomForests)

- Shepherd, J. C., et al. Emergence and Persistence of Cxcr4-Tropic Hiv-1 in a Population of Men from the Multicenter Aids Cohort Study. [J Infect Dis, 198, 1104-1112 (2008)](https://www.ncbi.nlm.nih.gov/pubmed/18783316).


==============
### Acknowledgments
This work made use of the High Performance Computing Resource in the Core Facility for Advanced Research Computing at Case Western Reserve University. We are thankful to Ms. Janet Schollenberger, Senior Project Coordinator, CAMACS, as well as Dr. Jeremy J. Martinson, Sudhir Penugonda, Shehnaz K. Hussain, Jay H. Bream, and Priya Duggal, for providing us the data related to the samples analyzed in the present study. Data in this manuscript were collected by the Multicenter AIDS Cohort Study (MACS) at (http://www.statepi.jhsph.edu/macs/macs.html) with centers at Baltimore, Chicago, Los Angeles, Pittsburgh, and the Data Coordinating Center: The Johns Hopkins University Bloomberg School of Public Health. The MACS is funded primarily by the National Institute of Allergy and Infectious Diseases (NIAID), with additional co-funding from the National Cancer Institute (NCI), the National Heart, Lung, and Blood Institute (NHLBI), and the National Institute on Deafness and Communication Disorders (NIDCD). MACS data collection is also supported by Johns Hopkins University CTSA. This study was supported by two grants from the National Institute of Health: NIDCR P01DE019759 (Aaron Weinberg, Peter Zimmerman, Richard J. Jurevic, Mark Chance) and NCI R01CA163739 (Hemant Ishwaran). The work was also partly supported by the National Science Foundation grant DMS 1148991 (Hemant Ishwaran) and the Center for AIDS Research grant P30AI036219 (Mark Chance). The contents of this publication are solely the responsibility of the authors and do not represent the official views of the granting agencies and institutions. The funders had also no role in the study design, data collection and analysis, decision to publish, or preparation of the manuscript.
