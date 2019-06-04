import os, sys

firstRun=input("first run number: ")
lastRun=input("Last run number: ")
runDuration=input("Run duration: ")

for iRun in range (firstRun, lastRun) :
    os.system('python acquire_sipm_laser_ext_trigger --config ../config/config.ini -o /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run%.4d --time %d' % (iRun, runDuration) )
