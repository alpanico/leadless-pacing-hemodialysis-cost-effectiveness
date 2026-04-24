# leadless-pacing-hemodialysis-cost-effectiveness
R code for the nationwide French medico-economic evaluation comparing leadless vs transvenous single-chamber pacemakers in haemodialysis patients. Analyses used linked SNDS/REIN data, propensity score weighting, and doubly robust net-benefit regression for survival and device-related infection outcomes.

# Cost-effectiveness of leadless versus transvenous cardiac pacing in haemodialysis patients

This repository contains the R code used for the analyses reported in the article:

**“Cost-effectiveness of leadless versus transvenous cardiac pacing in haemodialysis patients using nationwide administrative data.”**

## Purpose of this repository

The code is shared to document the analytical workflow used in the study, including:
- data management steps
- cost aggregation
- propensity score weighting
- doubly robust net-benefit regression analyses
- sensitivity analyses
- generation of tables and figures

## Data availability

The original individual-level data used in this study come from French nationwide administrative and registry sources (SNDS and REIN) and are subject to legal and regulatory restrictions.

These data cannot be made publicly available.

As a result, this repository is intended for **transparency of the analytical workflow**, not for full public reproducibility with the original data.

## Main analytical methods

The study compared leadless and transvenous single-chamber pacemakers in haemodialysis patients using:
- propensity score weighting
- doubly robust net-benefit regression
- simple-weighted and partitioned-time estimators
- overall survival and device-related infection-free survival as effectiveness outcomes

## Repository structure

- `code/`: R scripts used for the analyses
- `data/`: documentation only, no individual-level data
- `results/`: optional derived outputs that do not contain confidential information

## Software

Analyses were performed in **R**.

## Contact

Alexandre Panico : alexandre.panico@univ-lorraine.fr
