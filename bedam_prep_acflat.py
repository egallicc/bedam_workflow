# bedam python class
__doc__="""
$Revision: 0.1 $

A class to prepare BEDAM jobs using academic IMPACT

"""
# Contributors: Ahmet Mentes and Junchao Xia

import os, sys, time, re, glob
from operator import itemgetter, attrgetter
from schrodinger import structure, structureutil
from schrodinger.application import inputconfig
from schrodinger.utils import cmdline
import schrodinger.utils.log
import shutil
import signal
import glob
import time

from math import *

import sqlite3 as lite
from bedam_prep_ac import bedam_prep_ac


class bedam_prep_acflat (bedam_prep_ac):
    """
    Class to set up remd calculations with flattening 
    """
    def __init__(self, command_file, options):
        bedam_prep_ac.__init__(self, command_file, options)
#
# Impact input file templates for academic
#
    def resetTemplatesFlat(self):
        """ Setup templates for input files for academic impact"""

        self.input_remd = """
write file -
"%s" -
      title -
"%s" *

CREATE
  build primary name species1 type auto read sqldb file -
"%s"
  build primary name species2 type auto read sqldb file -
"%s"
QUIT

SETMODEL
  setpotential
    mmechanics nb12softcore umax %s %s consolv agbnp2
    weight constraints buffer %s
    weight bind rxid 0 nrep %d %s
%s
  quit
  read parm file -
"paramstd.dat" -
  noprint
  energy rest domain cmdist kdist %s dist0 %s toldist %s -
      read file "%s"
  energy parm dielectric 1 nodist -
   listupdate 10 -
    cutoff 12 hmass 5
  energy rescutoff byatom all
  zonecons auto
  energy constraints bonds hydrogens
  energy flatten torsions read file dihedral.dat
QUIT

MINI
  read restart coordinates formatted file "%s"
QUIT

DYNAMICS
  input cntl nstep %s delt 0.001
  input cntl constant temperature langevin relax 10.0 -
      rxmd nwalkers %d nstepexch %s -
      hamiltonian
  input target temperature %s
  input cntl initialize temperature at %s
  input cntl nprnt %s
  input cntl tol 1.00000e-07
  input cntl statistics off
  run rrespa fast 4
QUIT

DYNAMICS
  input cntl nstep %s delt 0.001
  input cntl constant temperature langevin relax 1.0
  input cntl rxmd nwalkers %d nstepexch 200 -
      hamiltonian %s
  input target temperature %s
  input cntl initialize temperature at %s
  input cntl nprnt %s
  input cntl tol 1.00000e-07
  input cntl statistics off
  write trajectory coordinates rxid every %s -
      external file "%s"
  run rrespa fast 4
  write restart coordinates rxid formatted file "%s"  
  write sql name species1 file "%s"
  write sql name species2 file "%s"
QUIT

END
"""

        self.input_remd_restart = """
write file -
"%s" -
      title -
"%s" *

CREATE
  build primary name species1 type auto read sqldb file -
"%s"
  build primary name species2 type auto read sqldb file -
"%s"
QUIT

SETMODEL
  setpotential
    mmechanics nb12softcore umax %s %s consolv agbnp2
    weight constraints buffer %s
    weight bind rxid 0 nrep %d %s
%s
  quit
  read parm file -
"paramstd.dat" -
  noprint
  energy rest domain cmdist kdist %s dist0 %s toldist %s -
      read file "%s"
  energy parm dielectric 1 nodist -
   listupdate 10 -
    cutoff 12 hmass 5
  energy rescutoff byatom all
  zonecons auto
  energy constraints bonds hydrogens
  energy flatten torsions read file dihedral.dat
QUIT

DYNAMICS
  input cntl nstep %s delt 0.001
  input cntl constant temperature langevin relax 1.0
  input cntl rxmd nwalkers %d nstepexch 200 -
      hamiltonian %s
  input target temperature %s
  input cntl initialize temperature at %s
  input cntl nprnt %s
  input cntl tol 1.00000e-07
  input cntl statistics off
  read restart coordinates rxid formatted file "%s"
  write trajectory coordinates rxid every %s -
      external file "%s"
  run rrespa fast 4
  write restart coordinates rxid formatted file "%s"  
  write sql name species1 file "%s"
  write sql name species2 file "%s"
QUIT

END
"""

    
##################### MAIN CODE ##########################
if __name__ == '__main__':

    # Setup the logger
    logger = schrodinger.utils.log.get_output_logger("bedam_prep_ac")

    # Parse arguments:
    usage = "%prog [options] <inputfile>"
    parser = cmdline.SingleDashOptionParser(usage)
    (options, args) = parser.parse_args(sys.argv[1:])
    
    if len(args) != 1:
        parser.error("Please specify ONE input file")
    
    commandFile = args[0]

    print ""
    print "========================================================="
    print "       BEDAM Job Preparation Using Academic IMPACT       "
    print "========================================================="
    print ""
    print "SCHRODINGER: " + os.environ['SCHRODINGER']
    print "Started at: " + str(time.asctime())
    print "Input file:", commandFile
    print ""
    sys.stdout.flush()
    
    print "Reading options"
    bedam = bedam_prep_acflat(commandFile, options)
    print "Reading the templates for input files ..."
    bedam.setupTemplates()
    print "Reading the reset templates flat for input files ..."
    bedam.resetTemplatesFlat()
    print "Get DMS files from maegz files ..."
    bedam.getDesmondDMSFiles()
    print "Writing BEDAM restraint file ..."
    bedam.writeRestraintFile()
    print "Adding atomic restraints to receptor file ..."
    bedam.writeRecStructureFile()
    print "Writing Minimization/Thermalization input file ..."
    bedam.writeThermInputFile()
    print "Writing Remd Production input file ..."
    bedam.writeRemdInputFile()
    print "Writing Trajectory conversion input file ..."
    bedam.writeReadTrajInputFile()

    print
    print "Job preparation complete"

#EOF
