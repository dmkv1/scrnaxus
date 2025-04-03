#!/bin/bash

SCRIPTDIR=$(dirname "$0")
SCRIPTPATH=$(realpath $SCRIPTDIR)

cd $SCRIPTPATH
docker build -t rstudio-server-scrnaseq:1.2 $SCRIPTPATH/rstudio-server-scrnaseq
docker compose up -d
