
****************************************************************************
This document provides detailed instructions on how to run the analysis 
code for the Jpsi polarization analysis using 2012 data via dielectron 
channel.


If you run into any problems with regards to using these code, please
contact Siwei Luo (sluo21@uic.edu)
****************************************************************************

****************************************************************************

PART I : StPicoDstMaker LBL PicoDst
The LBL PicoDst is created for general physics analysis that are being carried out by the RNC soft physics group. It is a ROOT-tree-like format data generated from the STAR standard MuDst data. The PicoDst data are produced, stored and archieved by this group.
All current source codes, scripts are maintained in the RCF CVS area. Here is a README file.
$CVSROOT/offline/users/dongx/pico/README
The latest source code location is:
$CVSROOT/offline/users/dongx/pico/source
All the source codes used for production for various data sets are located here:
$CVSROOT/offline/users/dongx/pico/prod
For more details on StPicoDstMaker, please refer to http://rnc.lbl.gov/~xdong/SoftHadron/picoDst.html.

------------------------------------------------------------------

PART II : Pico to AnaTree 
code directory is at PDSF /global/homes/h/huangbc/jpsi/pp200run12/pico2tree_evtAct
analysis.cxx, analysis_MB.cxx, analysis_sys.cxx, produce_Tree.cxx

They are copied to RCF /star/u/siwei/polarization/run12pol_review/PDSF_pico2tree_evtAct.

The output files are copied to RCF 
/star/u/siwei/polarization/run12pol_review/PDSF_pico2tree_evtAct/out_produce_1
and 
/star/u/siwei/polarization/run12pol_review/PDSF_pico2tree_evtAct/out_produce_2

a) The compile file is GNUmakefile

b) Execute the following commands to test the code and submit batch jobs
-cd /global/homes/h/huangbc/jpsi/pp200run12/pico2tree_evtAct/
- make
- To run a local test: ./analysis 
- To submit batch jobs for HT trigger: ./subjob_analysis_ht.sh
- To submit batch jobs for MB trigger: ./subjob_analysis_mb.sh

------------------------------------------------------------------

PART III : Read AnaTree and fill histogram of polarization of Jpsi in dataset
code directory is at PDSF /global/homes/s/siwei/link/test20160707newtree_2eid_inv30
analysis.cxx, analysis_MB.cxx

They are copied to RCF /star/u/siwei/polarization/run12pol_review/PDSF_newtree_2eid_inv30.

a) The compile file is GNUmakefile

b) Execute the following commands to test the code and submit batch jobs
-cd /global/homes/s/siwei/link/test20160707newtree_2eid_inv30
-make
-To run a local test: ./analysis 
-To submit batch jobs: ./subjobs.sh 

------------------------------------------------------------------

PART IV : Calculate efficiency of detectors

The codes are located in "StRoot/StMyElectronMaker" and "StRoot/StMyJpsiEffMakerSmearing66" : 
-StMyElectronMaker : the main analysis class for StMyElectron object
-StMyJpsiEffMakerSmearing66 : the analysis class for read tree and calculate the angles with regards to J/psi polarization.
a) Compile the codes
-cd /star/u/siwei/polarization/run12pol_review
-cons

b) Execute the following commands to test the code and submit batch jobs
-To run a local test: root -b -q doTest.C 
-To submit batch jobs: ./submit_jobs.sh

------------------------------------------------------------------

PART V: Toy Monte Carlo study
test3D_2eid_inv30_combine.C
a) This macro is used to do the Monte Carlo study of parameter estimator. The main goal is to check the bias of the estimator in both the parameter estimation and parameter statistical error estimation. The result is used to calibrate the estimator performance on real data fit. The value of parameters are feed to probability density function of J/psi decay products to generate the distribution of cos(theta) and phi. For each cos(theta) and phi bin in the 2D distribution histogram a random number is rolled and checked with the value in 2D efficiency histogram to decide whether fill the count into histogram to take the efficiency of detector into consideration and make the raw yield of cos(theta) phi distribution more like the real dataset. And then use the estimator to extract the paramters of probability density function of J/psi decay products. For each kinematic range of pT bin, five value of parameters with step 0.2 is measured from the estimator. A line is fit to the five measured value of parameters and parameters statistical error to correct the bias of estimator.  

b) Execute the macro
-cd /star/u/siwei/polarization/run12pol_review
-To run a local test: root -b -q test3D_2eid_inv30_combine.C
-To submit batch jobs: ./submit_test3D_2eid_inv30_jobs.sh 

------------------------------------------------------------------

PART VI: Fit the polarization parameters
sPlot3D_2eid_inv30_combine.C
a) The macro reads the dataset histogram on 2D cos(theta) phi distribution and 2D efficiency as a function of cos(theta) phi from embedding. Then calculate the log-likelihood and search for the minimum point of log-likelihood. And based on the toy MonteCarlo study of bias of estimator, correct both the central value estimation and statistical error of paramters accordingly. This macro is used to fit the negative log-likelihood function and extract polarization parameters' value. TMinuit package is applied for the parameter estimation.

b) Execute the macro to extract parameter
-cd test3D_2eid_inv30_combine.C
-To run a local test: root -b -q sPlot3D_2eid_inv30_combine.C
-To submit batch jobs: ./submit_sPlot3D_2eid_inv30_jobs.sh
This shell script is to organize the most macros used in the analysis including Monte Carlo study, the map between the measured parameters' value and true parameters' value, the invariant mass spectra, raw (cos(theta),phi) distribution , (cos(theta),phi) efficiency, the fit of J/psi polarization parameters, the calculation of systematic uncertainties and plotting the physics results of J/psi polarization parameters measurements. 

c) systematic uncertainty
sPlot3D_uncertainty.C
This macro is used to calculate the systematic uncertainty from different sources. The macro will read the output files on the parameters extraction. And calculate the maximum variation corresponding to different uncertainty source.  

------------------------------------------------------------------

PART VII : Plot the physics results

1. PWG_compare_invariant.C
a) This macro is used to plot the invariant mass spectra from dataset.
b) Execute the macro: root -b -q PWG_compare_invariant.C

2. PWGproposal.C
a) This macro is used to show the invariant mass spectra in different kinematic range and both the 2 dimensional raw (cos(theta),phi) distribution and 2 dimensional detector efficiency as a function of (cos(theta),phi). And the comparison between expected costheta and phi distribution and raw yield costheta and phi distribution from dataset.  
b) Execute the macro: root -b -q PWGproposal.C

3. PWG_Draw_Likelihood.C 
a) This macro is used to plot the minimum point position and the 1sigma contour of negative log-likelihood for one pT bin in a specific frame.
b) Execute the macro: root -b -q PWG_Draw_Likelihood.C 

4. PWG_All_Likelihood.C
a) This macro is used to plot the minimum point position and the 1sigma contour of negamtive log-likelihood for all pT bins and all frames.
b) Execute the macro: root -b -q PWG_All_Likelihood.C
c) submit_PWG_Draw_Likelihood.sh
This script is used to organize macros PWG_Draw_Likelihood.C and PWG_All_Likelihood.C excution. 

5. PWGMCvsData.C
a) This macro plots the J/psi polarization parameters and their statistical errors compared with the distribution of their counterparties result from Monte Carlo study. 
b) Execute the macro: root -b -q PWGMCvsData.C

6. uncertainty_pt.C
a) This macro is used to plot the physics results of J/psi polarization measurements including the J/psi polarization parameters' value, statistical error, systematic uncertainty, systematic uncertainty from different sources, parameters' value and statistical error before and after the estimator calibration, and the trend fit of J/psi polarization parameters.b) Execute the macro: root -b -q uncertainty_pt.C
b) Execute the macro: root -b -q uncertainty_pt.C




