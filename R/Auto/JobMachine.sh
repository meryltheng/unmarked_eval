#!/bin/bash

for i in {4..100}
do
  Rscript --vanilla ./R/Auto/2_generateDetections_Auto.R large $i
done