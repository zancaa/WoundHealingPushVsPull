#!/bin/bash
#
# This script assumes that the following has been run successfully:
# scons cl=0 co=1 b=GccOpt ts=projects/WoundHealingProliferativeHubLeadingEdge/test/TestBoundedGrowthModel.hpp
#

# 10 seeds
for ((i=1; i<=10; i++))
do

    echo "Beginning run ${i}."

    nice ../build/optimised/TestBoundedGrowthModelRunner -factive_min 0 -factive_max 10 -mqvf_min 0.6 -mqvf_max 1.0 -dmax_min 1 -dmax_max 21 -factive_num_sweeps 10 -mqvf_num_sweeps 8 -dmax_num_sweeps 20 -seed ${i} > output/bounded_${i}_output.txt 2>&1 &

done

echo "Jobs submitted"
