# Characterizing HIV-1 viral kinetics as a pharmacodynamic marker using nonlinear mixed-effects modeling 

A population pharmacodynamic analysis of HIV-1 viral load decline after the start of antiretroviral therapy, 
using a biphasic nonlinear mixed-effects (NLME) model in R.

## Overview

Understanding how quickly and how variably antiretroviral therapy suppresses HIV-1 is important for dose selection and regimen design. 
This project fits a biphasic NLME model to 361 HIV-1 RNA measurements from 46 patients in ACTG 315 (zidovudine, lamivudine and ritonavir), 
with baseline CD4 count as a covariate on the second-phase rate. The model is evaluated with diagnostic plots and a visual predictive check (VPC).

The data contain no drug concentrations, so this is a **pharmacodynamic** model. 
It does not support dosing conclusions, but it provides a structure to which pharmacokinetics and exposure-response relationships can be added.

## Model

Viral load (log10 copies/mL) follows a biphasic decay:

```
DV(t) = log10( exp[p1 - exp(ll1) * t] + exp[p2 - l2 * t] ) + e
```

- `p1`, `p2`: log-scale intercepts of the first and second phases
- `exp(ll1)`: first-phase decay rate (per day)
- `l2`: second-phase decay rate (per day), with baseline CD4 (centered, per 100 cells) as a covariate
- Independent random effects on all four parameters; additive residual error on the log10 scale
- Fitted by maximum likelihood with `nlme`

Values at or below the quantification limit (100 copies/mL, 11% of observations) are set to half the limit. The script also offers an `"exclude"` option for sensitivity analysis.

## Main results

| Parameter | Estimate (SE) | 
|---|---|
| p1 | 11.68 (0.20) | 
| ll1 | -1.19 (0.06) | 
| p2 | 6.34 (0.25) |
| l2, intercept | -0.0005 (0.0032)
| l2, CD4 effect | 0.0105 (0.0034) | 
| Between-patient SD | p1 1.11; ll1 0.26; p2 1.28; l2 0.017
| Residual SD | 0.38 | 


The VPC showed the observed 5th, 50th and 95th percentiles within their simulated 90% intervals. 

<img width="1608" height="795" alt="Screenshot 2026-10-22 at 2 41 57 PM" src="https://github.com/user-attachments/assets/1362d393-ca09-4d22-809e-383f1cfe40cd" />

See the full report for details and limitations.
[ACTG315_PD_Model_Report.pdf](https://github.com/user-attachments/files/33840938/ACTG315_PD_Model_Report.pdf)


## Repository structure

```
.
├── README.md
├── R/
│   └── viral_dynamics_analysis.R   # data preparation, model fit, diagnostics, VPC
├── data/                           # see Data
├── report/
│   └── ACTG315_PD_Model_Report.pdf # full report
```

## Data

The ACTG 315 data are available in this repository. Or download them from the dataset page maintained by Dr. Hulin Wu at UTHealth Houston School of Public Health:

https://sph.uth.edu/dept/bads/faculty-home/hulinwu/datasets/actg315longitudinaldatanlme

Save the file as `data/ACTG315.csv`. The script expects the columns `Patid`, `Day`, `log10.RNA.` and `CD4`. 
Please follow the data provider's terms of use.

## How to run

Requirements: R (4.x), and the packages `nlme` (3.1.168 was used), `dplyr` and `ggplot2`.

```r
install.packages(c("nlme", "dplyr", "ggplot2"))
source("R/viral_dynamics_analysis.R")
```

Settings you can change are at the top of the script: the quantification limit, the last day analyzed (196),
the censoring method (`"half_LLOQ"` or `"exclude"`), and the number of VPC simulations.

## Limitations

- No drug concentrations, so no PK/PD or dosing conclusions.
- Below-limit values are substituted not with censored method, which may bias second-phase estimates.
- The second-phase intercept is weakly identified, and all mechanistic interpretations are hypotheses.

## References

- Perelson AS, Essunger P, Cao Y, et al. Decay characteristics of HIV-1-infected compartments during combination therapy. *Nature*. 1997;387(6629):188-191. doi:10.1038/387188a0
- Lederman MM, Connick E, Landay A, et al. Immunologic responses associated with 12 weeks of combination antiretroviral therapy consisting of zidovudine, lamivudine and ritonavir: results of AIDS Clinical Trials Group Protocol 315. *J Infect Dis*. 1998;178:70-79.
- Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS*. Springer; 2000.
