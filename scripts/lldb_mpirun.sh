#!/bin/bash

if [ ${OMPI_COMM_WORLD_RANK} -eq "0" ]; then
    echo "Runnnig Debugger on node `hostname`"
    lldb $*
else
    $*
fi
exit 0
