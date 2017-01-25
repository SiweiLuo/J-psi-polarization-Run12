#!/bin/bash

date
#scheme="mle"

#read -p "1. rawdata; 2. efficiency; 3. getdata; 4. getefficiency; 5. update; 6. combine; 7. lambdas: " step
#read -p "choose scheme mle or x2: " scheme
#read -p "Rebin to 10*10 bins for likelihood histograms: " rebin
rebin=1

read -p "filemin : " filemin
read -p "filemax : " filemax

#file=0
#while [ $file -lt 100 ] 
#do
	ssh -YC siwei@rftpexp.rhic.bnl.gov "cd /star/u/siwei/polarization/test20160707/MLEmethod/ && ./generate_con_file.sh $filemin $filemax"	
	echo "Submit 2Dfit jobs on PDSF at $(date)"
	jobs_status="calculating"
	until [ "$jobs_status" == "DONE" ] 
	do
		jobs_status=$(ssh -YC siwei@rftpexp.rhic.bnl.gov 'set str=`condor_q -submitter siwei`
				if ("$str" == "") then
					echo "DONE"
				endif')
	done

	ssh -YC siwei@rftpexp.rhic.bnl.gov "cd /star/u/siwei/polarization/test20160707/MLEmethod/ && ./combine.sh $filemin $filemax && ./submit_fitlambda_jobs.sh $filemin $filemax"
#	ssh -YC siwei@pdsf.nersc.gov "cd /global/homes/s/siwei/polarization/test20160707/MLEmethod/ && ./submit_fitlambda_jobs.sh"
	echo "Submit 2DPOL8 jobs on PDSF at $(data)"
	jobs_status="calculating"
	until [ "$jobs_status" == "DONE" ] 
	do
		jobs_status=$(ssh -YC siwei@rftpexp.rhic.bnl.gov 'set str=`condor_q -submitter siwei`
				if ("$str" == "") then
					echo "DONE"
				endif')
	done
#	file=$(($file+1))
#done

echo " Calculation is finished at $(date) "
