#! /bin/sh

# This file is supposed to be ran from super4pcs_install_dir/scripts/
#
# run both 4pcs and super4pcs on the hippo model.
# parameters are optimised to get an accurate alignement of the two models
# -o input estimated overlap in [0,1]
# -d delta, used to compute the LCP between the two models
# -t maximum computation time in seconds
# -n number of samples used for the matching

export LD_LIBRARY_PATH=../lib/

time -p ../bin/Super4PCS -i ../assets/hippo1.obj ../assets/hippo2.obj -o 0.7 -d 0.01 -t 1000 -n 200 -r super4pcs_fast.obj -m mat_super4pcs_fast.txt
time -p ../bin/Super4PCS -i ../assets/hippo1.obj ../assets/hippo2.obj -o 0.7 -d 0.01 -t 1000 -n 200 -r      4pcs_fast.obj -m      mat_4pcs_fast.txt -x

