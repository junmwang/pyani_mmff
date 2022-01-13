#!/bin/csh
# Usage: ani.py  -i input.ani -o output.out -m method_id >run.log 
# method_id can be:
# 1-CG-BS,
# 2-CG-WS,
# 3-LBFGS,
# 4-LBFGS-WS,
# 5-BFGS,
# 0-single point
# In this tutorial, we have a mol2 file of alanine dipeptide with phi-psi angles at 210 and 270 degrees, 
# We want to constrained minimization with ani-2x

# Step 1: generate a gcrt file with antechamber from AMBERTools
antechamber -fi mol2 -fo gcrt -i ala2_210_270.mol2 -o ala2_210_270.gcrt -j 0

# Step 2: generate input file for running ani.py code. Note that the torsional angles defined by
# 5-7-9-15 and 7-9-15-17 will be frozen during minimization.
gcrt2ani -i ala2_210_270.gcrt -o ala2_210_270.ani -ct "5-7-9-15,7-9-15-17"

# Step 3: run ani.py to do minimization
$HOME/anaconda3/bin/python ani.py -i ala2_210_270.ani -o ala2_210_270.out -m 2 > run_ala2_210_270.log

# Step 4: extract minimized structure using ani2ani. 
ani2ani -i ala2_210_270.out -o ala2_210_270.pdb -e ala2_210_270.ene

# Step 5: this step is optional. In case a minimization is not converged, one may generate a restart file 
# using ani2ani
ani2ani -i ala2_210_270.out -c ala2_210_270.ani -o ala2_210_270_new.ani -f 2



