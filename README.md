		        BEDAM WORKFLOW  (v1.0)
		    =============================

			  Emilio Gallicchio
		      emilio.gallicchio@gmail.com

 Please acknowledge use of this software citing the following publication:

  Gallicchio E, Lapelosa M, Levy RM. Binding energy distribution
  analysis method (BEDAM) for estimation of protein-ligand binding
  affinities. J Chem Theory Comput, 6, 2961-2977 (2010).

 Copyright and Disclaimers
 -------------------------

 Copyright (c) 2010-2013 Emilio Gallicchio, emilio.gallicchio@gmail.com

 License Agreement and Disclaimers
 ---------------------------------

 This software is licensed under GPL v.3.
 
 http://www.gnu.org/licenses/gpl.html

 See LICENSE in this directory.

 Installation
 ------------

 Unpack the distribution in a directory of your choice. Dependencies
 includes numpy, and the MBAR python package of Michael Shirts and
 John Chodera downaloadable from http://simtk.org. These libraries
 should be installed where python can find them (for example by
 setting the PYTHONPATH env variable as appropriate)

 Introduction
 ------------

 The BEDAM Binding Energy Distribution Analysis Method is an absolute
 binding free energy estimation and analysis methodology based on a
 statistical mechanics theory of molecular association and efficient
 computational strategies built upon parallel Hamiltonian replica
 exchange, implicit solvation and multi-state statistical
 inference. The method takes its name from the technique it employs to
 extract standard binding free energies from the statistical analysis
 of the probability distributions of the energies of association over
 a series of conformational ensembles connecting the bound and unbound
 states. The ability to carry out extensive conformational sampling is
 one of the main advantages of BEDAM over existing FEP and absolute
 binding free energies protocols in explicit solvent which suffer from
 limited exploration of conformational space. The method has been
 extensively used to estimate protein-ligand and host-guest binding
 free energies.

 This python workflow facilitates the preparation and the analysis of
 BEDAM binding free energy calculations. It is designed to work with
 the IMPACT program within the Schrodinger computational
 environment. The workflow works in three steps:

1. System Preparation
 $SCHRODINGER/run bedam_prep.py bedam.cntl

2. Launch equilibration and productions calculations
$SCHRODINGER/impact -i <jobname>_mintherm.inp -LOCAL
$SCHRODINGER/impact -i <jobname>_remd.inp -LOCAL

3. Analysis
 $SCHRODINGER/run bedam_analyze.py bedam.cntl

A control file ('bedam.cntl' above) is used to specify parameters and
settings. Examples for some of the main settings are:

\#name of receptor .mae file
RECEPTOR_FILE 'bcy_noprop.maegz'
\#name of ligand .mae file
LIGAND_FILE 'benzene.maegz'
\#list of lambdas for each replica (16 replicas, in this case). 
\#lambda=0 is the decoupled state, lambda=1 is the coupled state.
\#The binding free energy is the free energy difference between the lambda=1
\#and lambda=0 states plus a standard state correction term.
LAMBDAS '0.0,0.001,0.002,0.004,0.005,0.006,0.008,0.01,0.02,0.04,0.07,0.1,0.25,0.5,0.75,1.0'
\#the atoms of the receptor and ligand that define their centroids. 
\#These are given in ASL (atom selection language). A flat-bottom
\#harmonic restraint (see below) on the centroid-to-centroid distance defines
\#the binding site region and keeps the ligand and the receptor together at
\#small lambdas.
REST_LIGAND_CMRECASL '( all) AND NOT (( atom.ele H ) )'
REST_LIGAND_CMLIGASL '( all) AND NOT (( atom.ele H ) )'
\#parameters of the flat-bottom harmonic restraints between the centroids
\#defined above.
\# force constant in kcal/mol/A^2 of receptor-ligand restraints
REST_LIGAND_CMKF 3.0
\# equilibrium distance in Angstroms of receptor-ligand restraint
REST_LIGAND_CMDIST0 0.0
\# distance tolerance in Angstroms of receptor-ligand restraint
REST_LIGAND_CMTOL 6.0
\#the atoms of the receptor that are harmonically restrained
REST_RECEPTOR_ASL '(not (atom.num 3, 5, 10, 11, 14, 16, 21, 22, 25, 27, 32, 33, 36, 38, 43, 44, 47, 49, 54, 55, 58, 60,
65, 66, 69, 71, 76, 77 ) ) AND NOT (( atom.ele H ) )'
\# force constant of receptor atomic restraints in kcal/mol/A^2
REST_RECEPTOR_KF 0.6
\#Temperature in K
TEMPERATURE 300
\# number of equilibration MD steps
EQUILIBRATION_STEPS 10000
\# number of production MD steps.
PRODUCTION_STEPS 500000
\#Frequency of printing information in output file in number of steps. The number
\#of binding energy samples collected for each replica is
\#PRODUCTION_STEPS/PRNT_FREQUENCY = 500 in this case.
\#500 x 16replicas = 800 samples in total
PRNT_FREQUENCY 1000
\# Frequency of saving trajectory frames. 500 frames per replica.
TRJ_FREQUENCY 1000

The workflow takes as input .mae files of receptor and ligand, plus a
definition (see above) of the binding site region in terms of
a receptor-ligand restraint potential. Analysis of the results produce,
among other things, the estimated values of the binding free energy.

See the 'examples' directory for an example of usage for a host guest
system. Similar procedures are used for protein-ligand receptors.
