PMA-ES: Probabilistic Morphological Analysis with Entropy Screening
Replication material for:
> Faye, B., Cardebat, J.-M. *A coherent map of the possible: entropy screening and copula-calibrated dependence for probabilistic foresight.* Submitted to *Technological Forecasting and Social Change*.
Corresponding author: Benoît Faye (bfaye@inseec.com), INSEEC Grande École Bordeaux.
Contents
File	Description
`PMA-ES_reproduction.R`	Reproduces all results of the paper: the real case (Tables 3 to 6, Janus factors, sensitivity to the sign threshold and to ρmax) and the simulation study (Table 2, Figures 1 to 3).
`PMA-ES_analyst_tool.R`	Applies the method to any probabilistic morphological study through the function `generate_families()`.
`results_PMA-ES.xlsx`	Full results produced by the reproduction script, one sheet per table or analysis.
Requirements
R (version 4.0 or later). Both scripts run in base R with no external dependency. The packages `writexl` or `openxlsx` are optional and only used to export results to Excel. Without them, results are written as CSV files.
Reproducing the results
```r
source("PMA-ES_reproduction.R")
```
The full run takes about 5 to 10 minutes, most of it for the simulation study (500 replications per condition). Lower `REPS_SIM` at the top of the script for a quick test. Random seeds are fixed, so the run reproduces the figures reported in the paper. Copula probabilities are estimated by Monte Carlo with 10⁶ draws.
Applying the method to your own study
```r
source("PMA-ES_analyst_tool.R")
res <- generate_families(experts = my_experts, mdi = my_mdi)
```
The method needs two materials that a probabilistic morphological study already produces: the direct influence matrix (MDI, rated 0 to 3) and the probabilities each expert assigns to the hypotheses of each factor. The header of the script describes the input format, all parameters with their default values, and the outputs. A demonstration on the real case at the end of the script shows the exact data format.
Data
The real case reanalyzes the foresight study published in Faye, B., Cardebat, J.-M., Hocquard, H. (2023). Shared stakeholder views on the future of the world wine industry competitiveness. Systèmes alimentaires / Food Systems, 8, 137–165. Expert responses are anonymized.
Use
The code and data are made available for research purposes. Please cite the paper above if you use them.
