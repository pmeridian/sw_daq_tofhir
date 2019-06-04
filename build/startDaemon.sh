rm -v /tmp/d.sock
rm -v /dev/shm/daqd_shm
./daqd --socket-name=/tmp/d.sock --daq-type=GBE
