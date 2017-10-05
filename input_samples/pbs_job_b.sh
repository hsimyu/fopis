#!/bin/sh
#QSUB -q {QUEUE}
#QSUB -ug {GROUP}
#QSUB -W 24:00
#QSUB -A p=1:t=1:c=1
#QSUB -m e

bin_name=tdpic
path_to_bin=../${bin_name}

if [ -e ${path_to_bin} ]; then
    mpiexec.hydra ${path_to_bin}  > ./std.out
    RC=$?
    if [ $RC -gt 0 ]; then
        echo "[ERROR] return code: $RC"
    fi
    exit $RC
else
    echo "[ERROR] An executionable ${bin_name} file does not exist at ${path_to_bin}."
    exit
fi
