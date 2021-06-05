FROM rocker/r-ver:3.6.1

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  libxml2-dev \
  libcairo2-dev \
  libsqlite-dev \
  libmariadbd-dev \
  libmariadbclient-dev \
  libpq-dev \
  libssh2-1-dev \
  unixodbc-dev \
  libsasl2-dev \
 ## from rocker/rstudio:3.6.1 Dockerfile
  file \
  git \
  libapparmor1 \
  libcurl4-openssl-dev \
  libedit2 \
  libssl-dev \
  lsb-release \
  psmisc \
  procps \
  python-setuptools \
  sudo \
  wget \
  libclang-dev \
  libclang-3.8-dev \
  libobjc-6-dev \
  libclang1-3.8 \
  libclang-common-3.8-dev \
  libllvm3.8 \
  libobjc4 \
  libgc1c2 \
  libudunits2-dev \
  libgdal-dev \
  libxt-dev \
  && install2.r --error \
    --deps TRUE \
    tidyverse \
    dplyr \
    devtools \
    formatR \
    remotes \
    selectr \
    caTools \
    BiocManager \
    githubinstall \
    rstan \
    ggmcmc \
    brms