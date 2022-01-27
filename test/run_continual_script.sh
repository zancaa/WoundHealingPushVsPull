#!/bin/bash
# 
# This script assumes that the following has been run successfully:
# scons cl=0 co=1 b=GccOpt ts=projects/WoundHealingProliferativeHubLeadingEdge/test/TestContinualGrowthModel.hpp
#

# 10 seeds
for ((i=1; i<=10; i++))
do

    echo "Beginning run ${i}."

    nice ../build/optimised/TestContinualGrowthModelRunner -factive_min 0 -factive_max 10 -mqvf_min 0.6 -mqvf_max 1.2 -factive_num_sweeps 10 -mqvf_num_sweeps 12 -seed ${i} > output/continual_${i}_output.txt 2>&1 &

done

echo "Jobs submitted"
