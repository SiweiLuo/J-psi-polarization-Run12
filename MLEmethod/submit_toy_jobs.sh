#!/bin/bash
date 
rebin=$1

for file in {200..299}
do
logdir=~/polresults/20160707/MLEmethod/log/sys_${file}
if [ ! -d $logdir ] ;then 
	mkdir $logdir
fi

rundir=~/polresults/20160707/MLEmethod/log/run_${file}
if [ ! -d $rundir ] ;then 
	mkdir $rundir
fi

comdir=~/polresults/20160707/MLEmethod/combinedrootfiles/sys_${file}
if [ ! -d $comdir ] ;then 
	mkdir $comdir
fi

mledir=~/polresults/20160707/MLEmethod/rootfiles/sys_${file}
if [ ! -d $mledir ] ; then 
	mkdir $mledir
#	for trig in {0..2}
#	do
#		for xsect in 1 11 21 31 41 51 61 71 81 91 
#		do
#			for ysect in 1 11 21 31 41 51 61 71 81 91 
#			do
#				echo "$date" >> $logdir/toy${file}_${trig}_${xsect}_${ysect}.log 
#			   	echo "$date" >> $logdir/toy${file}_${trig}_${xsect}_${ysect}.err
#				qsub -hard -l h_vmem=2.0G -l projectio=1 -o $logdir/toy${file}_${trig}_${xsect}_${ysect}.log -e $logdir/toy${file}_${trig}_${xsect}_${ysect}.err ./mlemethod.sh $file $trig $xsect $ysect $rebin

#				condor_submit  

#				echo "Job# $jobid submitted ..."
#			done
#		done
#	done
fi
done 
