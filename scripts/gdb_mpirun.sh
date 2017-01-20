#!/bin/bash

if [ ${OMPI_COMM_WORLD_RANK} -eq "0" ]; then
    echo "Runnnig GDB on node `hostname`"
    xterm -e gdb --args $*
else
    $*
fi
exit 0
