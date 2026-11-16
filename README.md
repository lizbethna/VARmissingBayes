##--------------------------------------------------

Paper: A Bayesian vector autoregressive model to address missing data and prediction in precision beekeeping

Authors: M. Carmen Robustillo (1), Lizbeth Naranjo (2), M. Isabel Parra (3), Carlos J. Perez (3)

(1) Universidad de Málaga, Spain. 

(2) Departamento de Matematicas, Facultad de Ciencias, Universidad Nacional Autonoma de Mexico, Ciudad de Mexico, Mexico. 

(3) Departamento de Matematicas, Universidad de Extremadura, Caceres, Spain.

Journal: Applied Soft Computing Journal

Year: 2026

##--------------------------------------------------

ABSTRACT:
Precision beekeeping increasingly depends on sensor-based monitoring, yet these systems frequently generate incomplete multivariate time series due to connectivity interruptions or temporary hardware failures. These gaps highlight the importance of using statistical models that manage missing data, ensuring accurate predictions and timely interventions. This work proposes a Bayesian vector autoregressive model with exogenous variables that integrates imputation and forecasting in a unified probabilistic framework. Missing observations are treated as latent variables and estimated jointly with model parameters via Hamiltonian Monte Carlo, ensuring coherent uncertainty propagation and preserving the temporal structure of hive dynamics. A simulation study found that model parameter estimates remain stable across varying levels of missing data, with no significant differences compared to estimates from the no-missing-data scenario. When applied to real data from the BeeObserver project, the proposed approach provided comparable accuracy rates when imputing data with different proportions of missing values. Additionally, it provides reliable predictions for key beekeeping variables, achieving a mean absolute error in one (seven)-day forecast of less than 0.85 (2.3)$^\circ C$ for temperature, 3.2 (4)\% for humidity, and 420 (710) g for weight under a cross-validation procedure. The proposed approach thus provides a coherent probabilistic basis for integrating forecasting and data reconstruction in practical hive monitoring systems, supporting informed decision-making even in the presence of extended sensor outages. 

Keywords:
Bayesian statistics; Hamiltonian Monte Carlo; Missing data; Precision beekeeping; Sensor data; Time series prediction; Vector autoregressive model.

##--------------------------------------------------

Instructions to run the codes in Markdown, R and Stan software. 
The codes are applied to obtain one similar analysis with cross-validations as in Section:
5. Experimental results

The code is implemented for a rolling-window cross-validation procedure, trained on a 90-days window and to produce forecast for subsequent 7 days, getting 11 windows. 
PDF is run for one window, but R code can be modified for run the complte cross-validation.

##--------------------------------------------------
Files: 

Bee135 data (this is not the original dataset):
- Bee135.txt

VAR(p) model: 
- Bee135_BVARX2_thw.Rmd;
- Bee135_BVARX2_thw.pdf; 
- BVARX_p.stan


VAR(p)-missing data model:
- Bee135_BVARX2_missing_thw.Rmd;
- Bee135_BVARX2_missing_thw.pdf; 
- BVARX_p_missing.stan 

##--------------------------------------------------

To run the files, do the following.
 
1.- Install the packages necessary to run the R files. These are indicated in the R file. 

2.- Change the address indicated in 'address="HERE"' and ‘setwd("HERE")’. 
This is the address where the data and Stan files are in.

3.- Run the code in file "Bee135_BVARX2_xxx.Rmd" to estimate the parameters of the  model.

##--------------------------------------------------

This research was funded by Ministry of Science, Innovation and Universities - Spain, and State Research Agency - Spain 
(Projects PID2021-122209OB-C32 y PID2024-155179NB-C21), 
funded by MICIU/AEI/10.13039/501100011033 and European Union (European Regional Development Fund).


[LogoSolo.pdf](https://github.com/user-attachments/files/23808555/LogoSolo.pdf)

##-------------------------------------------------- 
