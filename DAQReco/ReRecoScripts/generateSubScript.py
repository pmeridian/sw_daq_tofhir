#! /usr/bin/env python

import os
import sys
import optparse
import datetime
import time

#JOB PARAMS
EXEC        = 'runReco.sh'
BINARY      = 'ConvertTOFPETSinglesToEvents'
EOSFOLDERNAME = '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RecoData/v2/RecoWithTracks/'
RUNRANGE    = []
NJOBS       = 1

f = open("/eos/home-m/mtd/www/BTL/MTDTB_FNAL_Feb2020/triggerTimeCalibration/fitParameters.txt","r")
lines = f.readlines()
for word in lines:
    if "#" in word:
        continue
    RUNRANGE.append(word.split()[0])
f.close()

pwd = os.environ['PWD']
EXEC = pwd+"/"+EXEC
BINARY = pwd+"/"+BINARY
current_time = datetime.datetime.now()
timeMarker = "submit_%04d%02d%02d_%02d%02d%02d" % (current_time.year,current_time.month,current_time.day,current_time.hour,current_time.minute,current_time.second)
workingDir = pwd+"/"+timeMarker
os.system("mkdir -p "+workingDir)


#prepare condor sub fileName
with open(workingDir+"/condor.sub", "w") as fo:
    fo.write("+JobFlavour = \"espresso\"\n\n")
    fo.write("transfer_output_files = \"\"\n")
    fo.write("executable = "+EXEC+"\n")
    fo.write("log        = "+workingDir+"/output.log\n")
    fo.write("requirements = (OpSysAndVer =?= \"CentOS7\")\n")
    fo.write("max_retries = 3\n")
    fo.write("\n")

    for run in RUNRANGE:
        fo.write("output = "+workingDir+"/"+run+".out\n")
        fo.write("error = "+workingDir+"/"+run+".err\n")
        fo.write("arguments = "+run+" "+EOSFOLDERNAME+" "+BINARY+"\n")
        fo.write("queue "+str(NJOBS)+"\n")
        fo.write("\n")
