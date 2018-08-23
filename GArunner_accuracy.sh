#!/bin/bash
chmod -R 755 pyGAsubmit_accuracy.sh
for graphfilename in *.gpickle; do
	chmod -R 755 $graphfilename;
	sbatch pyGAsubmit_accuracy.sh $graphfilename;
done