# Overview

MATLAB, R, and RStan code for omics-based metabolic flux estimation without labeling for extended trans-omic analysis (OMELET).

OMELET is an approach to use simultaneously obtained multi-omic data to infer metabolic fluxes in each condition, identify changes in metabolic flux between conditions, and quantify contributions of regulators to the changes in metabolic flux. 



## Workflow

The series of analysis performed in Uematsu et al. are described in `run_OMELETmouse.m` for mouse data and `run_OMELETyeast.m` for simulated data from the yeast kinetic model.

1. Make input data for OMELET (MATLAB).
   + `make_input_OMELET.m` for mouse data.
   + `@simulateKineticModel` for simulated data from the yeast kinetic model.
2. Perform OMELET to infer metabolic fluxes and other parameters (RStan).
   + `run_OMELETmouse.R` for mouse data.
   + `run_OMELETyeast.R` for simulated data from the yeast kinetic model.
3. Calculate contributions of regulators to changes in metabolic flux between conditions (MATLAB).
   + `prep` method in `@outputOMELETmouse` class for mouse data.
4. Make figures (MATLAB and R).
   + `makeFig*` methods in `@outputOMELETmouse` class for mouse data.
   + Figs. 3A and S3 are plotted by R. 

![](OMELET_graphical_model.png)

# Requirements

The developmental version of the programs has been tested on the following systems:

Ubuntu 18.04.5 LTS

MATLAB R2019a

SimBiology toolbox (for application of OMELET to the yeast kinetic model)

Docker version 20.10.6

# Contents

## 1. Make input data for OMELET

### make_input_OMELET

MATLAB scripts to make input data for OMELET.

+ `make_input_OMELETmouse.m`

  This code makes input data for OMELET based on `S_OMELETmouse.csv`.

+ `S_OMELETmouse.csv`

  Stoichiometric matrix and information on cofactors and allosteric effectors for `OMELETmouse.stan` model.

+ `S_OMELETyeast.csv`

  Stoichiometric matrix and information on cofactors and allosteric effectors for `OMELETyeast.stan` model.

### @simKineticModel

A MATLAB class to generate datasets (metabolites, enzymes, and metabolic fluxes) in several conditions from the yeast kinetic model as input data for OMELET. To make input data for OMELET, this class load `S_OMELETyeast.csv`.



## 2. Perform OMELET to infer metabolic fluxes and other parameters

### OMELET_rstan

Scripts to perform OMELET in RStan.

+ `run_OMELETmouse.R`

  This code performs parameter estimation of `OMELETmouse.stan` model. You need to prepare input data `input_OMELETmouse` by `make_input_OMELETmouse.m`. This code also makes Figure 3A.

+ `run_OMELETyeast.R`

  This code performs parameter estimation of `OMELETyeast.stan` model. You need to prepare input data in `input_OMELETyeast` by `@simKineticModel/makeRstanInput.m`. This code also makes Figure S3.



The RStan environment can be build from docker image `saori/rstan` from [DockerHub](https://hub.docker.com/r/saori/rstan). Make sure you have [Docker Engine](https://docs.docker.com/engine/install/) installed.

+ `saori/rstan:latest`

  This docker image is to run RStudio Server (based on `rocker/rstudio:3.6.1`). Dockerfile is provided in `docker_files_for_rstan/saori_rstan_latest/Dockerfile`. 

  To get `saori/rstan:latest`:

  ```shell
  docker pull saori/rstan:latest
  ```

  To run this image in a docker container:

   ```shell 
   bash run_docker_RStan.sh container_name rstudio
   ```

  

+ `saori/rstan:cmd`

  This docker image is to run R console without running RStudio Server (based on `rocker/r-ver:3.6.1`). You can use this image to perform parameter estimation of several models in one computer. Dockerfile is provided in `docker_files_for_rstan/saori_rstan_cmd/Dockerfile`.

  To get `saori/rstan:cmd`:
  
  ```shell
  docker pull saori/rstan:cmd
  ```
  
  To run this image in a docker container:
  
  ```shell
  bash run_docker_RStan.sh container_name cmd
  ```
  



## 3. Calculate contributions of regulators to changes in metabolic flux between conditions

## 4. Make figures

### @outputOMELETmouse

A MATLAB class to make outputs of OMELET. This calculates contributions of regulators to changes in metabolic flux between conditions, as well as make figures for Uematsu et al. You need to prepare RStan input by `make_input_OMELETmouse.m` and output by `run_OMELETmouse.R`.



## data

Data used for Uematsu et al.

+ `metabolome.csv`, `proteome.csv`, `transcriptome.csv`

  Metabolomic, proteomic, and transcriptomic data from wild-type (WT) and leptin-deficient obese (*ob*/*ob*) mice at fasting state and after oral glucose administration. These are used as input data for OMELET (`make_input_OMELETmouse.m`) and to make figures (`@outputOMELETmouse`).

+ `blood_glucose.csv`, `blood_insulin.csv`

  Time-course data of blood glucose and insulin after oral glucose administration. These data are used only for make Figure S1 in Uematsu et al.

+ `pathway_matlab1.csv`

  A template pathway of the glucose metabolism for constructing trans-omic networks.

+ `BIOMD0000000503_url.xml`

  The yeast kinetic model downloaded from [BioModels](https://www.ebi.ac.uk/biomodels/BIOMD0000000503) for generating datasets (metabolites, enzymes, and metabolic fluxes) in several conditions as input data for OMELET (`@simulateKineticModel`).



## +lib

Shared program libraries for the above scripts.



# Contact

Saori Uematsu: suematsu@bs.s.u-tokyo

Satoshi Ohno: sohno@bs.s.y-tokyo.ac.jp

Shinya Kuroda: skuroda@bs.s.u-tokyo.ac.jp

# Reference

Saori Uematsu, Satoshi Ohno, Kaori Tanaka, Atsushi Hatano, Toshiya Kokaji, Yuki Ito, Hiroyuki Kubota, Ken-ichi Hironaka, Yutaka Suzuki, Masaki Matsumoto, Keiichi I. Nakayama, Akiyoshi Hirayama, Tomoyoshi Soga, and Shinya Kuroda. Omics-based label-free metabolic flux inference reveals dysregulation of glucose metabolism in liver associated with obesity.

