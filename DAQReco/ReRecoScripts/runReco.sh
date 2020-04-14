#!/bin/sh

RUN=${1}
OUTDIR=${2}
EXEC=${3}

NAME=run${RUN}_events.root
WORKDIR=`pwd`

echo "Run number is ${RUN}"
echo "Work dir: ${WORKDIR} at `hostname`"

mkdir -p ${OUTDIR}

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh

${EXEC} ${RUN}
mv outFile_${RUN}.root ${OUTDIR}/${NAME}
