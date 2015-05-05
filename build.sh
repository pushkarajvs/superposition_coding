#!/bin/bash

echo "Building the superposition coding system"
make
echo "Make successful"

echo "Now running the executables"

echo "First, generating near user 12 rates"
./single_user_12rates_near_op > output_near.txt
echo "Far ran perfectly"

echo "Next, generating far user 12 rates"
./single_user_12rates_far_op > output_far.txt
echo "Near ran perfectly"

echo "Finally, finding the extra rate-pairs feasible due to superposition coding"
./superposition_coding_op > output_spc.txt
echo "Generated the extra rate-pairs"

echo "Plotting all the rates and also the convex hull obtained by the superposition coding"
octave plot_rate_region.m
echo "Finished successfully"
