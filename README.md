# Overview

MATLAB, R, and RStan code for omics-based metabolic flux estimation without labeling for extended trans-omic analysis (OMELET).

OMELET is an approach to use simultaneously obtained multi-omic data to infer metabolic fluxes in each condition, identify changes in metabolic flux between conditions, and quantify contributions of regulators to the changes in metabolic flux. 

![](OMELET_graphical_model.png)


# Requirements

The developmental version of the programs has been tested on the following systems:

Ubuntu 18.04.5 LTS or Windows 10

MATLAB R2019a

SimBiology toolbox (for application of OMELET to the yeast kinetic model)

Docker version 20.10.6


## Docker environment

The RStan environment can be build from docker images `saori/rstan` from [DockerHub](https://hub.docker.com/r/saori/rstan). Make sure you have [Docker Engine](https://docs.docker.com/engine/install/) installed. We provide two options.

+ `saori/rstan:latest`

  This docker image is used to run RStudio Server (based on `rocker/rstudio:3.6.1`). Dockerfile is provided as `docker_files_for_rstan/saori_rstan_latest/Dockerfile`. 

  To get `saori/rstan:latest`:

  ```shell
  docker pull saori/rstan:latest
  ```


+ `saori/rstan:cmd`

  This docker image is used to run R console without running RStudio Server (based on `rocker/r-ver:3.6.1`). You can use this image to perform parameter estimation of several models in one computer. Dockerfile is provided as `docker_files_for_rstan/saori_rstan_cmd/Dockerfile`.

  To get `saori/rstan:cmd`:
  
  ```shell
  docker pull saori/rstan:cmd
  ```
  
  To run this image in a docker container:
  
  (terminal in Ubuntu)
  
  ```shell
  cd OMELET_rstan
  bash run_docker_RStan.sh container_name cmd
  ```
  
  (command prompt in Windows):
  
  ```
  cd OMELET_rstan
  run_docker_RStan.bat container_name cmd
  ```
  


# Application to mouse data

## Step 1. Make input data for OMELET

First, you need to make input data for OMELET. Execute the following commands on MATLAB;

To use the functions under `./make_input_OMELET`:
```matlab
addpath <path to ./OMELET/make_input_OMELET>
```

To specify the path to the .csv file that contains stoichiometric matrix and information on cofactors and allosteric effectors in a specific format.
```matlab
S_path = <path to ./OMELET/make_input_OMELET/S_OMELETmouse.csv>:
```

To specify the path to the directory containing multi-omic data (metabolome.csv, proteome.csv, transcriptome.csv):
```matlab
data_dir_path = <path to ./OMELET/data>`
```

To specify the directory to save the input data for OMELET:
```matlab
savedir = <path to the directory where you want to save your input data>;
```

To choose indexes from `Index` row in metabolome.csv, proteome.csv, transcriptome.csv:
```matlab
idx_smplgrp = [1 2 3 4];
% 1: WT in the fasting state
% 2: WT 4h after oral glucose administration
% 3: ob/ob in the fasting state
% 4: ob/ob 4h after oral glucose administration
```

After specifying all the arguments, execute the following commands on MATLAB.

```matlab
make_input_OMELETmouse(S_path,data_dir_path,idx_smplgrp,savedir);
```

The input data for OMELET will be saved under `savedir`, and those used in Uematsu et al. are under `./OMELET/OMELET_rstan/input_OMELETmouse_sample`.


## Step 2. Perform OMELET to infer metabolic fluxes and other parameters

### Prepare your working directory

Create your working directory in your computer.

Copy the following functions under `./OMELET/OMELET_rstan` to your working directory:
```
initf_OMELETmouse.R
my_plot.R
my_rstan_OMELETmouse.R
my_rstan_opt.R
OMELETmouse.stan
run_docker_RStan.sh
```

Create a directory for input data, and copy all the input data generated in step 1 under this directory.


### Run a Docker container

`cd` to your working directory.

To run the `saori/rstan:latest` image in a docker container;

+ If you are using Ubuntu, execute the following commands on the terminal:
```shell 
bash run_docker_RStan.sh <your container name> rstudio
```

+ If you are using Windows, execute the following commands on the command prompt (not on the PowerShell):
```
run_docker_RStan.bat <your container name> rstudio
```

After running your Docker container, visit `localhost:8787` in your browser and log in with username `rstudio` and password `1`. The Rstudio environment will be launched on your browser where Rstan environments have already been installed.


### Run OMELET






To plot Figure 4A, execute the following command on R.
```R

```


## Step 3. Calculate contributions of regulators to changes in metabolic flux between conditions

After moving to the cloned main project directory `./OMELET`, execute the following commands on MATLAB;

To specify the path to the directory containing multi-omic data (metabolome.csv, proteome.csv, transcriptome.csv):
```matlab
data_dir_path = <path to ./OMELET/data>
```

To specify the path to the directory containing the results of parameter estimation of OMELET by Rstan.
```matlab
rstan_path = <path to the results by Rstan>
```

To specify the path to the .mat object that is generated by `make_input_OMELETmouse.m` in step 1.
```matlab
model_path = <path to the 'model_data.mat'>;
```

After specifying all the arguments, execute the following commands on MATLAB.
```matlab
obj = outputOMELETmouse(rstan_path,model_path,data_path);
obj.prep;
```
The resulting `obj` object contains definitions of metabolic network, omics data, inferred metabolic fluxes, the contributions of regulators to changes in metabolic flux between conditions. 


## Step 4. Make figures

First you need to specify the directory to save figures.
```matlab
savedir = <path to the directory where you want to save your figures>;
mkdir(savedir);
```

Using the `obj` object created in step 3, execute the following commands on MATLAB.
```matlab
obj.makeFig2(savedir);
obj.makeFig4B(savedir);
obj.makeFig4C_4D(savedir);
obj.makeFig5B_7B_S8B(savedir);
obj.makeFig5C_7A_S8A(savedir);
obj.makeFig6BCD_7C_S8C(savedir);
obj.makeFigS1(savedir);
obj.makeFigS4(savedir);
obj.makeFigS5(savedir);
obj.makeFigS6(savedir);
obj.makeFigS7(savedir);
```

To save the `obj` object under `savedir`, execute the following command.
```matlab
save([savedir '/obj_outputOMELETmouse.mat'],'obj');
```


# Application to yeast data

For the analysis of simulated data from the yeast kinetic model:

1. Make input data for OMELET (MATLAB).
   + `@simulateKineticModel` class
2. Perform OMELET to infer metabolic fluxes and other parameters (RStan).
   + `OMELET_rstan/run_OMELETyeast.R` 
   + Figure S3 is plotted by R in this script.


## data

Data used for Uematsu et al.

+ `data/metabolome.csv`, `data/proteome.csv`, `data/transcriptome.csv`

  Metabolomic, proteomic, and transcriptomic data from wild-type (WT) and leptin-deficient obese (*ob*/*ob*) mice at fasting state and after oral glucose administration. These are used as input data for OMELET (`make_input_OMELET/make_input_OMELETmouse.m`) and to make figures (`@outputOMELETmouse`).

+ `data/pathway_matlab1.csv`

  A template pathway of the glucose metabolism for constructing trans-omic networks.

+ `data/BIOMD0000000503_url.xml`

  The yeast kinetic model downloaded from [BioModels](https://www.ebi.ac.uk/biomodels/BIOMD0000000503) for generating datasets (metabolites, enzymes, and metabolic fluxes) in several conditions as input data for OMELET (`@simulateKineticModel`).



Additional data not used in OMELET but for Figure S1 in Uematsu et al.

+ `data/blood_glucose.csv`, `data/blood_insulin.csv`

  Time-course data of blood glucose and insulin after oral glucose administration.



## +lib

Shared program libraries for the above scripts.



# Contact

Saori Uematsu: suematsu@bs.s.u-tokyo.ac.jp, saorid6e9p6p9@gmail.com

Satoshi Ohno: sohno@bs.s.u-tokyo.ac.jp

Shinya Kuroda: skuroda@bs.s.u-tokyo.ac.jp

# Reference

Saori Uematsu, Satoshi Ohno, Kaori Tanaka, Atsushi Hatano, Toshiya Kokaji, Yuki Ito, Hiroyuki Kubota, Ken-ichi Hironaka, Yutaka Suzuki, Masaki Matsumoto, Keiichi I. Nakayama, Akiyoshi Hirayama, Tomoyoshi Soga, and Shinya Kuroda. Omics-based label-free metabolic flux inference reveals dysregulation of glucose metabolism in liver associated with obesity.

