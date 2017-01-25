#!/bin/bash

cp scripts/runMc.con scripts/rawMcAll.con

filemin=$1
filemax=$2
rebin=1
#for file in {0..19}
for file in $(seq ${filemin} ${filemax})
do
	for trig in {0..2}
	do
		for pt in {0..5}
		do
			for frame in {0..1}
			do
				#generate MC raw data
				cp scripts/run.csh scripts/raw_${file}_${trig}_${pt}_${frame}.csh
				echo "nohup root -b << EOF ">>scripts/raw_${file}_${trig}_${pt}_${frame}.csh
				echo ".x rawMCdata.C($file,$trig,$pt,$frame)">>scripts/raw_${file}_${trig}_${pt}_${frame}.csh
				echo ".q">>scripts/raw_${file}_${trig}_${pt}_${frame}.csh
				echo "EOF">>scripts/raw_${file}_${trig}_${pt}_${frame}.csh
	
				chmod 755 scripts/raw_${file}_${trig}_${pt}_${frame}.csh

				echo "Executable    = scripts/raw_${file}_${trig}_${pt}_${frame}.csh">>scripts/rawMcAll.con
				echo "Output        = /star/u/siwei/polresults/20160707/MLEmethod/log/run_${file}/raw_${file}_${trig}_${pt}_${frame}.out">>scripts/rawMcAll.con
				echo "Error         = /star/u/siwei/polresults/20160707/MLEmethod/log/run_${file}/raw_${file}_${trig}_${pt}_${frame}.err">>scripts/rawMcAll.con
				echo "Queue"        >> scripts/rawMcAll.con
				echo "   	       ">>scripts/rawMcAll.con
			
				echo ${file} ${trig} ${pt} ${frame}
				#generate MC raw data
			done
		done
	done 
done

#chmod 755 scripts/raw*.csh

condor_submit scripts/rawMcAll.con


