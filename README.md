# IRSF
Interaction Random Survival Forest (IRSF): Ensemble Survival Tree Models to Reveal Variable Interactions in Association with Time-to-Events Outcomes


==============
### License

PRIMsrc is open source / free software, licensed under the GNU General Public License version 3 (GPLv3), 
sponsored by the [Free Software Foundation](http://www.fsf.org/). To view a copy of this license, visit 
[GNU Free Documentation License](http://www.gnu.org/licenses/gpl-3.0.html).


==============
### Description
R codes for the analyses and the generation of the results shown in manuscript (Dazard et al., 2017). Codes contain randomization, interaction modeling and prediction subroutines to be used in addition to the following R packages (http://cran.r-project.org/): [`survival`](https://CRAN.R-project.org/package=survival) for Kaplan-Meier and Cox regression modeling, [`NADA`](https://CRAN.R-project.org/package=NADA) for correlation analysis in the presence of censoring, [`randomForestSRC`](https://CRAN.R-project.org/package=randomForestSRC) for RSF modeling (Ishwaran and Kogalur, 2013, 2007), and [`ggRandomForests`](https://CRAN.R-project.org/package=ggRandomForests) (Ehrlinger, 2014) for Random Forrest exploration/visualization. Default parameter specifications were used for all main functions.

==============
### References
- Dazard J-E., Ishwaran H., Mehlotra R.K., Weinberg A. and Zimmerman P.A. Ensemble Survival Tree Models to Reveal Variable Interactions in Association with Time-to-Events Outcomes. 
[submitted (2017)]

- Ishwaran, H. and Kogalur, U.B. Contributed R Package Random Forests for Survival, Regression and Classification (RF-SRC).
[CRAN (2013)](https://CRAN.R-project.org/package=randomForestSRC)

- Ishwaran H. and Kogalur U.B. Random Survival Forests for R. 
[R News, 7(2), 25-31, (2007)](https://pdfs.semanticscholar.org/951a/84f0176076fb6786fdf43320e8b27094dcfa.pdf)

- Ehrlinger J. Contributed R Package ggRandomForests: Visually Exploring Random Forests
[CRAN (2014)](https://CRAN.R-project.org/package=ggRandomForests)
