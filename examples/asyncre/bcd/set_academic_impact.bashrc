#!/bin/sh
IMPACTHOME=/home/emilio/schrod/academic-impact/
export IMPACTHOME

IMP_ROOT=$IMPACTHOME
export IMP_ROOT

IMPACT_EXEC=$IMP_ROOT/bin/Linux-x86_64
export IMPACT_EXEC

LD_LIBRARY_PATH=$IMP_ROOT/lib/Linux-x86_64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
