#!/bin/bash

# $1: container's name
# $2: whether launch rstudio server

if [ "$2" = "rstudio" ]; then
  CURRENT=$(cd $(dirname $0);pwd)
  docker run -it \
             -v $CURRENT:/home/rstudio \
             -p 8787:8787 \
             -e PASSWORD=1 \
             --name $1 \
             saori/rstan:latest
elif [ "$2" = "cmd" ]; then
  CURRENT=$(cd $(dirname $0);pwd)
  docker run -it \
             -v $CURRENT:/home/rstudio \
             -w /home/rstudio \
             --name $1 \
             saori/rstan:cmd
else
   echo "invalid argument: please specify rstudio or cmd."
fi
