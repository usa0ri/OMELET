@echo off

if "%~2"=="rstudio" (
   docker run -it -v "%CD%":/home/rstudio -p 8787:8787 -e PASSWORD=1 --name %1 saori/rstan:latest
) else if "%~2"=="cmd" (
   docker run -it -v "%CD%":/home/rstudio -w /home/rstudio --name %1 saori/rstan:cmd
) else (
  echo "invalid argument: please specify rstudio or cmd."
)