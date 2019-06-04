RUN=$1
echo "trigger:"
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 06 00" | wc | awk '{print "00 00 06 00: "$1}'
echo "single bars:"
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 00 00" | wc | awk '{print "00 00 00 00: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 00 01" | wc | awk '{print "00 00 00 01: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 00 06" | wc | awk '{print "00 00 00 06: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 00 08" | wc | awk '{print "00 00 00 08: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 00 14" | wc | awk '{print "00 00 00 14: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 00 15" | wc | awk '{print "00 00 00 15: "$1}'
echo "array:"
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 00" | wc | awk '{print "00 00 02 00: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 01" | wc | awk '{print "00 00 02 01: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 02" | wc | awk '{print "00 00 02 02: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 03" | wc | awk '{print "00 00 02 03: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 04" | wc | awk '{print "00 00 02 04: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 05" | wc | awk '{print "00 00 02 05: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 06" | wc | awk '{print "00 00 02 06: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 07" | wc | awk '{print "00 00 02 07: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 08" | wc | awk '{print "00 00 02 08: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 09" | wc | awk '{print "00 00 02 09: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 10" | wc | awk '{print "00 00 02 10: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 11" | wc | awk '{print "00 00 02 11: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 12" | wc | awk '{print "00 00 02 12: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 13" | wc | awk '{print "00 00 02 13: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 14" | wc | awk '{print "00 00 02 14: "$1}'
./print_raw -i /data/TOFHIR_MTDTB_FNAL_Jun2019/RAW/run$1 | grep "00 00 02 15" | wc | awk '{print "00 00 02 15: "$1}'

