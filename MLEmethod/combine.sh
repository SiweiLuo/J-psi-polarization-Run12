#!/bin/bash
filemin=$1
filemax=$2

for file in $(seq $filemin $filemax)
do
	intdir=~/polresults/20160707/MLEmethod/rootfiles/sys_${file}
	outdir=~/polresults/20160707/MLEmethod/combinedrootfiles/sys_${file}

	if [ ! -d $outdir ] ; then 
		mkdir $outdir
	fi

	rm $outdir/lambda_trg*.root

#	for trig in {0..2}
	for trig in 2
	do
#		for pt in {0..5}
		for pt in 3
		do
#			for frame in 0 1 
			for frame in 0 
			do
				hadd -f $outdir/lambda_trg${trig}_pt${pt}_frame${frame}.root $intdir/likelihood_${file}_${trig}_${pt}_${frame}_* 	
			done 
		done
	done
done

