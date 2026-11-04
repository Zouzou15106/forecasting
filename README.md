# PMA-ES: Probabilistic Morphological Analysis with Entropy Screening

Replication material 

> Contents

|File|Description|
|-|-|
|`PMA-ES\_reproduction.R`|Reproduces all results of the paper: the real case (Tables 4 to 7, Janus factors, sensitivity to the sign threshold, to ρmax and to the choice of copula) and the simulation study (Tables 2 and 3, Figures 1 to 3).|
|`PMA-ES\_analyst\_tool.R`|Applies the method to any probabilistic morphological study through the function `generate\_families()`.|
|`results\_PMA-ES.xlsx`|Full results produced by the reproduction script, one sheet per table or analysis.|

## Requirements

R (version 4.0 or later). Both scripts run in base R with no external dependency. The packages `writexl` or `openxlsx` are optional and only used to export results to Excel. Without them, results are written as CSV files.

## Reproducing the results

```r
source("PMA-ES\_reproduction.R")
```

The full run takes about 30 to 60 minutes, most of it for the two simulation studies (500 replications per condition). Lower `REPS\_SIM` and `REPS\_BENCH` for a quick test. The benchmark uses two cores on Linux and macOS and one on Windows. Each replication has its own seed, so the results do not depend on the number of cores. Random seeds are fixed, so the run reproduces the figures reported in the paper. Copula probabilities are estimated by Monte Carlo with 10⁶ draws.

## Applying the method to your own study

```r
source("PMA-ES\_analyst\_tool.R")
res <- generate\_families(experts = my\_experts, mdi = my\_mdi)
```

The method needs two materials that a probabilistic morphological study already produces: the direct influence matrix (MDI, rated 0 to 3) and the probabilities each expert assigns to the hypotheses of each factor. The header of the script describes the input format, all parameters with their default values, and the outputs. A demonstration on the real case at the end of the script shows the exact data format.

## Data

The real case reanalyzes the foresight study published in Faye, B., Cardebat, J.-M., Hocquard, H. (2023). Shared stakeholder views on the future of the world wine industry competitiveness. *Systèmes alimentaires / Food Systems*, 8, 137–165. Expert responses are anonymized.

## Use

The code and data are made available for research purposes. Please cite the paper above if you use them.

