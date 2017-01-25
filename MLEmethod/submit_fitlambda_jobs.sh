#!/bin/bash

date
#file=$1
filemin=$1
filemax=$2

cp scripts/runMc.con scripts/runFitAll.con

for file in $(seq ${filemin} ${filemax})
do
	logdir=~/polresults/20160707/MLEmethod/log/sys_${file}
	if [ -d $logdir ];then
		mkdir $logdir
	fi
	pdfdir=~/polresults/20160707/MLEmethod/pdf/sys_${file}
	if [ -d $pdfdir ];then
		mkdir $pdfdir
	fi

#	for trig in {0..2}
	for trig in 2
	do
#		for pt in {0..5}
		for pt in 3
		do
#			for frame in {0..1}
			for frame in 0
			do
#				qsub -hard -l h_vmem=2.0G -l projectio=1 -o $logdir/fitlambda${file}_${trig}_${pt}_${frame}.log -e $logdir/fitlambda${file}_${trig}_${pt}_${frame}.err ./fitlambda.sh $file $trig $pt $frame
#				echo "Job# $jobid submitted ..."
				cp scripts/run.csh scripts/fitlambda_${file}_${trig}_${pt}_${frame}.csh	
				
				echo "nohup root -b << EOF ">>scripts/fitlambda_${file}_${trig}_${pt}_${frame}.csh
				echo ".x fitlambda.C($file,$trig,$pt,$frame)">>scripts/fitlambda_${file}_${trig}_${pt}_${frame}.csh
	            echo ".q">>scripts/fitlambda_${file}_${trig}_${pt}_${frame}.csh
	            echo "EOF">>scripts/fitlambda_${file}_${trig}_${pt}_${frame}.csh

				echo "Executable    = scripts/fitlambda_${file}_${trig}_${pt}_${frame}.csh">>scripts/runFitAll.con
				echo "Output       	= /star/u/siwei/polresults/20160707/MLEmethod/log/run_${file}/fitlambda_${file}_${trig}_${pt}_${frame}.out">>scripts/runFitAll.con
				echo "Error         = /star/u/siwei/polresults/20160707/MLEmethod/log/run_${file}/fitlambda_${file}_${trig}_${pt}_${frame}.err">>scripts/runFitAll.con
#				echo "Log   = /star/u/siwei/polresults/20160707/MLEmethod/log/run_${file}/fitlambda_${file}_${trig}_${pt}_${frame}.log">>scripts/runFitAll.con
				echo "Queue" >>scripts/runFitAll.con
				echo "         ">>scripts/runFitAll.con
			done
		done
	done
done

chmod 755 scripts/fitlambda*.csh

condor_submit scripts/runFitAll.con



