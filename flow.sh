#!/bin/bash

date
read -p "1. rawdata; 2. efficiency; 3. getdata; 4. getefficiency; 5. update; 6. combine; 7. lambdas: " step

if [ "$step" = "step1" ]; then
	read -p "Recalculate the raw data? [Y/N]" yn 
	case $yn in 
		y ) echo "Recalculate the raw data at $(date)";; 
		N ) echo "Raw data remains the same";;
	esac  
fi

if [ "$step" = "step1" ] || [ "$step" = "step2" ];then
	read -p "Recalculate the efficiency? [Y/N]" yn
	case $yn in 
		y ) echo " Recalculate the efficiency at $(date)";;
		N ) echo " Efficiency remains the same";;
	esac  
fi

if [ "$step" = "step1" ] || [ "$step" = "step2" ] || [ "$step" = "step3" ]; then
	echo " Download raw data and efficiency rootfiles at $(date) "
	read -p "Download raw data? [Y/N]" yn
	case $yn in 
		y )
			strdate=$(date +%Y%m%d)
			if [ ! -d rootfile_${strdate} ]; then
				mkdir rootfile_${strdate}
			fi
			ssh -YC siwei@pdsf.nersc.gov '/global/homes/s/siwei/link/test20160707newtree/output/./addhistogram.sh'

			if [ ! -d rootfile_back_${strdate} ]; then 	
				mkdir rootfile_back_${strdate}
			fi
	
			cp rootfile/*.root rootfile_backup_${strdate}
			echo "Rootfiles are backed up at:
				$(date)" >> rootfile_backup_${strdate}/ReadMe.txt 
			echo "Rootfiles are backed up at:
				$(date)" >> rootfile/ReadMe.txt
	
			scp siwei@pdsf.nersc.gov:/global/homes/s/siwei/link/test20160707newtree/output/out_result1_20160919_ht*_1TrkPid_*/ht*_trg*_sys?.ana.root rootfile/ && scp siwei@pdsf.nersc.gov:/global/homes/s/siwei/link/test20160707newtree/output/out_result1_20160919_ht*_1TrkPid_*/ht*_trg*_sys??.ana.root rootfile/

			scp 'siwei@rftpexp.rhic.bnl.gov:/star/data01/pwg/siwei/Jpsi/test20160706/rootfile/OutFile_sys*.root' rootfile/ 

			echo "Rootfiles are updated at:
				$(date)" >> rootfile/ReadMe.txt ;;
	esac
fi

if [ "$step" = "step1" ] || [ "$step" = "step2" ] || [ "$step" = "step3" ] || [ "$step" = "step4" ];then
	echo " Scan lambdas on PDSF at $(date) " 
	read -p "Update the polarization directory on PDSF? [Y/N]" yn
	case $yn in 
		y ) rsync -au ~/polarization/ siwei@pdsf.nersc.gov:~/polarization/	;;
	esac 
fi

if [ "$step" = "step1" ] || [ "$step" = "step2" ] || [ "$step" = "step3" ] || [ "$step" = "step4" ] || [ "$step" = "step5" ];then
file=0
while [ $file -lt 1 ] 
do
	ssh -YC siwei@pdsf.nersc.gov "cd /global/homes/s/siwei/polarization/test20160707/ && ./submit_pol_jobs.sh $file"	
	echo "Submit jobs on PDSF at $(date)"
	jobs_status="calculating"
	until [ "$jobs_status" == "DONE" ] 
	do
		jobs_status=$(ssh -YC siwei@pdsf.nersc.gov 'set str=`qstat -u siwei`
				if ("$str" == "") then
					echo "DONE"
				endif')
	done
	ssh -YC siwei@pdsf.nersc.gov "cd /global/homes/s/siwei/polarization/test20160707/ && ./combine.sh $file"
	file=$(($file+1))
done
echo "Jobs are DONE at $(date)"
fi



if [ "$step" = "step1" ] || [ "$step" = "step2" ] || [ "$step" = "step3" ] || [ "$step" = "step4" ] || [ "$step" = "step5" ] || [ "$step" = "step6" ];then
echo " Download results from PDSF at $(date) " 

	rsync -au siwei@pdsf.nersc.gov:~/polresults/ ~/polresultspdsf/ 

	root -b  -q lambdas.C 
fi
echo " Calculation is finished at $(date) "
