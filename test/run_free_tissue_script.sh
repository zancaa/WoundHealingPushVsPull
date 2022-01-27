#!/bin/bash
# 
# This script assumes that the following has been run successfully:
# scons cl=0 co=1 b=GccOpt ts=projects/WoundHealingProliferativeHubLeadingEdge/test/TestFreeTissue.hpp
#

echo "Beginning simulation."

nice ../build/optimised/TestFreeTissueRunner -factive_min 0 -factive_max 10 -factive_num_sweeps 10 > output/free_tissue_output.txt 2>&1 &
