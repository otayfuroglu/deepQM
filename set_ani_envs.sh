#! /usr/bin/env bash

export PATH=/usr/local/cuda-9.2/bin${PATH:+:${PATH}}
export LD_LIBRARY_PATH=/usr/local/cuda-9.2/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

###### NC_ROOT EXPORTS ########
export NC_ROOT="/truba/home/otayfuroglu/opt/ASE_ANI"
export LD_LIBRARY_PATH="$NC_ROOT/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$NC_ROOT/lib:$PYTHONPATH"

