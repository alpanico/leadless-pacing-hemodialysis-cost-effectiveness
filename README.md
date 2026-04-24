# Cost-effectiveness of leadless versus transvenous cardiac pacing in haemodialysis patients

This repository contains the R code used for the analyses reported in the article:

**“Cost-effectiveness of leadless versus transvenous cardiac pacing in haemodialysis patients using nationwide administrative data.”**

## Overview

This project reports the analytical workflow used for a nationwide French medico-economic evaluation comparing leadless and transvenous single-chamber pacemakers in haemodialysis patients.

The analyses were based on linked French administrative and registry data (SNDS/REIN) and used propensity score weighting together with doubly robust net-benefit regression to evaluate cost-effectiveness for:
- overall survival
- device-related infection-free survival

## Purpose of this repository

The code is shared to document the main analytical workflow used in the study, including:
- propensity score construction and weighting
- weighted cost description
- doubly robust net-benefit regression models
- sensitivity analyses
- generation of tables and figures

Because the original data are confidential, this repository is intended for **transparency of the analytical approach**, not for full public reproducibility using the original individual-level data.

## Data availability

The original individual-level data used in this study come from French nationwide administrative and registry sources (SNDS and REIN) and are subject to legal and regulatory restrictions.

These data cannot be made publicly available.

The scripts provided here assume access to secured analytic datasets that were created in the authorized environment.

## Main analytical methods

The study compared leadless and transvenous single-chamber pacemakers in haemodialysis patients using:
- propensity score weighting
- doubly robust net-benefit regression
- simple-weighted and partitioned-time estimators
- multiple cost perspectives
- overall survival and device-related infection-free survival outcomes

## Repository structure

- `code/`: R scripts used for the analyses
- `data/`: documentation only; no individual-level data are shared
- `results/`: optional non-confidential derived outputs

## Software

Analyses were performed in **R**.

## Contact

Alexandre Panico  
alexandre.panico@univ-lorraine.fr
