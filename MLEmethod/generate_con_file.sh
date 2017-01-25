#!/bin/bash
date
cp scripts/runMc.con scripts/runMcAll.con

filemin=$1
filemax=$2
rebin=1
witheff=1
#for file in {0..19}
for file in $(seq ${filemin} ${filemax})
do
#	for trig in {0..2}
	for trig in 2
	do
#		for pt in {0..5}
#		do
#			for frame in {0..1}
#			do
				#generate MC raw data
#				echo "Executable    = scripts/rawMCdata_${file}_${trig}_${pt}_${frame}.csh">>scripts/runMcAll.con
#				echo "Error         = /star/u/siwei/polresults/20160707/MLEmethod/log/run_${file}/raw_${file}_${trig}_${pt}_${frame}.err">>scripts/runMcAll.con
#				echo "Log           = /star/u/siwei/polresults/20160707/MLEmethod/log/run_${file}/raw_${file}_${trig}_${pt}_${frame}.log">>scripts/runMcAll.con
#				echo "Queue"        >> scripts/runMcAll.con
#				echo "        ">>scripts/runMcAll.con
				#generate MC raw data
#			done
#		done
		for xsect in 1 11 21 31 41 51 61 71 81 91 
		do
			for ysect in 1 11 21 31 41 51 61 71 81 91 
			do

				cp scripts/run.csh scripts/run_${file}_${trig}_${xsect}_${ysect}.csh

				echo "nohup root -b << EOF ">>scripts/run_${file}_${trig}_${xsect}_${ysect}.csh
				echo ".x mlemethod.C($file,$trig,$xsect,$ysect,$rebin,$witheff)">>scripts/run_${file}_${trig}_${xsect}_${ysect}.csh
				echo ".q">>scripts/run_${file}_${trig}_${xsect}_${ysect}.csh
				echo "EOF">>scripts/run_${file}_${trig}_${xsect}_${ysect}.csh

				chmod 755 scripts/run_${file}_${trig}_${xsect}_${ysect}.csh
#mv run_$file.csh tasklist/.

				echo "Executable    = scripts/run_${file}_${trig}_${xsect}_${ysect}.csh">>scripts/runMcAll.con
				echo "Output        = /star/u/siwei/polresults/20160707/MLEmethod/log/run_${file}/run_${file}_${trig}_${xsect}_${ysect}.out">>scripts/runMcAll.con
				echo "Error         = /star/u/siwei/polresults/20160707/MLEmethod/log/run_${file}/run_${file}_${trig}_${xsect}_${ysect}.err">>scripts/runMcAll.con
#				echo "Log           = /star/u/siwei/polresults/20160707/MLEmethod/log/run_${file}/run_${file}_${trig}_${xsect}_${ysect}.log">>scripts/runMcAll.con
				echo "Queue"        >> scripts/runMcAll.con
				echo "        ">>scripts/runMcAll.con
				
				echo ${file} ${trig} ${xsect} ${ysect}
#    let "file=file+1"

			done
		done
	done 
done

echo "scripting is done"
#chmod 755 scripts/run*.csh

echo "submitting jobs"
condor_submit scripts/runMcAll.con


