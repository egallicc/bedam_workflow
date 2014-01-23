# bedam python class
__doc__="""
$Revision: 0.1 $

A class to prepare and analyze BEDAM jobs

"""
# Contributors: Emilio Gallicchio

import os, sys, time, re, glob
from operator import itemgetter, attrgetter
from schrodinger import structure, structureutil
from schrodinger.application import inputconfig
from schrodinger.utils.fileutils import get_structure_file_format, splitext, is_poseviewer_file, get_next_filename
from schrodinger.utils import cmdline
import schrodinger.utils.log
import shutil
import signal
import glob

from math import *
from numpy import * # numerical array library
import MBAR


# Setup the logger
logger = schrodinger.utils.log.get_output_logger("bedam_prep")

class bedam_job:
    """
    Class to set up BEDAM calculations
    """
    def __init__(self, command_file, options):
        self.command_file = command_file
        self.jobname = os.path.splitext( os.path.basename(command_file) )[0]
        self.impact = os.path.join(os.environ['SCHRODINGER'],'impact')

        self.setupTemplates()

        self.parseInputFile()

        self.printStatus() #debug: write out the parameters, etc.

    def error(self, text):
        """ Print an error line to the log file """
        logger.error(text)
        sys.stdout.flush()
    
    def exit(self, text):
        """ Print an error and exit """
        self.error(text)
        print 'exiting...'
        sys.exit(1)

    def parseInputFile(self):
        """
        Read keywords from BEDAM input file
        """
        config = inputconfig.InputConfig(self.command_file)
        self.keywords = dict(config)

    def printStatus(self):
        print 'command_file=', self.command_file
        print 'jobname=', self.jobname
        for k, v in self.keywords.iteritems():
            print k, v

#
# Impact input file templates
#
    def setupTemplates(self):
        self.input_idx =  """
write file -
"%s" -
      title -
"%s" *

CREATE
  build primary name species1 type auto read maestro file -
"%s"
  build primary name species2 type auto read maestro file -
"%s"
  build types name species1
  build types name species2
QUIT

SETMODEL
  setpotential
    mmechanics consolv agbnp2
  quit
  read parm file -
"paramstd.dat" -
  noprint
  energy parm dielectric 1 nodist -
   listupdate 10 -
    cutoff 12 hmass 5
  energy rescutoff byatom all
  zonecons auto
  energy constraints bonds hydrogens
QUIT

MINIMIZE
  write maestro name species1 file "%s"
  write maestro name species2 file "%s"
QUIT

END
"""
        self.input_mintherm = """
write file -
"%s" -
      title -
"%s" *

CREATE
  build primary name species1 type auto read maestro file -
"%s"
  build primary name species2 type auto read maestro file -
"%s"
  build types name species1
  build types name species2
QUIT

SETMODEL
  setpotential
    mmechanics consolv agbnp2
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
QUIT

MINIMIZE
  conjugate dx0 5.000000e-02 dxm 1.000000e+00
  input cntl mxcyc 200 rmscut 1.000000e-02 deltae 1.000000e-07
  run
QUIT

put 100 into 'temp0'
put %s into 'tempt'
put 10 into 'n'
put 'tempt'- 'temp0' into 'dtemp'
put 'dtemp' / 'n' into 'dt'

put 0 into 'i'
while 'i' lt 'n'

DYNAMICS
  input cntl nstep 1000 delt 0.0005
  input cntl constant totalenergy
  input cntl initialize temperature at 'temp0'
  input cntl nprnt 100
  input cntl tol 1.00000e-07
  input cntl stop rotations
  input cntl statistics off
  run rrespa fast 8
QUIT

put 'temp0' + 'dt' into 'temp0'
put 'i' + 1 into 'i'

endwhile

DYNAMICS
  write restart coordinates formatted file "%s"  
  write maestro file "%s"
QUIT

END
"""
        self.input_remd = """
write file -
"%s" -
      title -
"%s" *

CREATE
  build primary name species1 type auto read maestro file -
"%s"
  build primary name species2 type auto read maestro file -
"%s"
  build types name species1
  build types name species2
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
  write maestro file -
   "%s"
QUIT

END
"""

        self.input_remd_restart = """
write file -
"%s" -
      title -
"%s" *

CREATE
  build primary name species1 type auto read maestro file -
"%s"
  build primary name species2 type auto read maestro file -
"%s"
  build types name species1
  build types name species2
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
  write maestro file -
   "%s"
QUIT

END
"""



        self.input_reservoir = """ -
    l0reservoir nreservoir %s resvrxid 0 -
    traj nfile 2 maxrec %s nskip 1 delt 0.001 -
    coordinates temperature every 1 -
    traj fnames external file "%s" -
                         file "%s" -
    cmrestraints """

  
        self.input_trj = """

write file -
"%s" -
      title -
"%s" *

CREATE
  build primary name species1 type auto read maestro file -
"%s"
  build primary name species2 type auto read maestro file -
"%s"
  build types name species1
  build types name species2
QUIT

SETMODEL
  setpotential
    mmechanics nb12softcore umax %s consolv agbnp2
    weight constraints buffer %s
    weight bind rxid 0 nrep %d %s
  quit
  read parm file -
"paramstd.dat" -
  noprint
  energy parm dielectric 1 nodist -
   listupdate 10 -
    cutoff 12 hmass 5
  energy rescutoff byatom all
  zonecons auto
  energy constraints bonds hydrogens
QUIT

put 0 into 'repl'
while 'repl' lt %d

put $q$ concat (char 'repl') into 'bs3'
put 'bs3' concat $.trj$ into 'trjname'

show 'trjname'

TABLE
traj nfile 1 maxrec 50000000 nskip 1 delt 0.001 -
coordinates rxid every 1 -
traj fnames external file -
'trjname'
starttrack
QUIT

reset 'current.rxid'

show 'repl'
show 'current.rxid'

put $q$ concat (char 'current.rxid') into 'ms3'
put 'ms3' concat $.maegz$ into 'maename'

MINIMIZE
   conjugate dx0 5.000000e-02 dxm 1.000000e+00
   input cntl mxcyc 0 rmscut 1.000000e-02 deltae 1.000000e-07
   run
   write maestro file 'maename' append
QUIT

TABLE
stoptrack
QUIT

put 'repl' + 1 into 'repl'
endwhile

END
"""

        self.input_trj_rescore = """

write file -
"%s" verbose 1 -
      title -
"%s" *

CREATE
  build primary name species1 type auto read maestro file -
"%s"
  build primary name species2 type auto read maestro file -
"%s"
  build types name species1
  build types name species2
QUIT

SETMODEL
  setpotential
    mmechanics nb12softcore umax %s consolv agbnp2
    weight constraints buffer %s
    weight bind rxid 0 nrep %d %s
  quit
  read parm file -
"paramstd.dat" -
  noprint
  energy parm dielectric 1 nodist -
   listupdate 10 -
    cutoff 12 hmass 5
  energy rescutoff byatom all
  zonecons auto
  energy constraints bonds hydrogens
QUIT

put 0 into 'repl'
while 'repl' lt %d

put $q$ concat (char 'repl') into 'bs3'
put 'bs3' concat $.trj$ into 'trjname'

show 'trjname'

TABLE
traj nfile 1 maxrec 50000000 nskip 1 delt 0.001 -
coordinates rxid every 1 -
traj fnames external file -
'trjname'
starttrack
QUIT

reset 'current.rxid'

show 'repl'
show 'current.rxid'

if 'current.rxid' eq %d
MINIMIZE
   conjugate dx0 5.000000e-02 dxm 1.000000e+00
   input cntl mxcyc 0 rmscut 1.000000e-02 deltae 1.000000e-07
   run
QUIT
endif

TABLE
stoptrack
QUIT

put 'repl' + 1 into 'repl'
endwhile

END
"""



#
# runs Impact on the input mae files to get versions with internal 
# atom indexes
#
    def getImpactMaeFiles(self):
        receptor_file =  self.keywords.get('RECEPTOR_FILE')
        if not receptor_file:
            msg = "bedam_prep: No receptor file specified in the input file"
            self.exit(msg)
        if not os.path.exists(receptor_file):
            msg = 'File does not exist: %s' % receptor_file
            self.exit(msg)
        ligand_file =  self.keywords.get('LIGAND_FILE')
        if not ligand_file:
            msg = "bedam_prep: No ligand file specified in the input file"
            self.exit(msg)
        if not os.path.exists(ligand_file):
            msg = 'File does not exist: %s' % ligand_file
            self.exit(msg)
        impact_input_file =   self.jobname + '_idx' + '.inp'
        impact_output_file =  self.jobname + '_idx' + '.out'
        impact_jobtitle =     self.jobname + '_idx'
        out_receptor_file =   self.jobname + '_rcpt_idx' + '.maegz'
        out_ligand_file =     self.jobname + '_lig_idx' + '.maegz'
        f = open(impact_input_file, 'w')
        input =  self.input_idx % (impact_output_file, impact_jobtitle, receptor_file, ligand_file, out_receptor_file, out_ligand_file )
        f.write(input)
        f.close()
        cmd = self.impact + " -i " + impact_input_file + " -WAIT"
        os.system(cmd)
        if not os.path.exists(out_receptor_file) or not os.path.exists(out_ligand_file):
            print " failed"
            msg = "Impact job to analyze structure files failed"
            self.exit(msg)
        self.recidxfile = out_receptor_file
        self.ligidxfile = out_ligand_file


#
#  write the binding site restraint file *_cmrestraint.dat 
#
    def writeRestraintFile(self):
        # check that structure files with internal indexes have been generated
        if self.recidxfile is None or self.ligidxfile is None:
            msg = "writeRestraintFile: Internal error: structure files not found"
            self.exit(msg)
             
        receptor_asl =  self.keywords.get('REST_LIGAND_CMRECASL')
        if not receptor_asl:
            msg = "bedam_prep: No ASL specified for receptor site center"
            self.exit(msg)
        ligand_asl =  self.keywords.get('REST_LIGAND_CMLIGASL')
        if not ligand_asl:
            ligand_asl = '( all)'

        try:
            st = structure.StructureReader(self.recidxfile).next()
        except:
            print "Cannot open Maestro structure file %s" % self.recidxfile
            sys.exit(-1)
        nrecatoms = len(st.atom)
        atoms = structureutil.evaluate_asl(st, receptor_asl)
        rec_atoms = [];
        for iat in atoms:
            rec_atoms.append(st.atom[iat].property['i_i_internal_atom_index'])

        #computes receptor CM
        cmrx = cmry = cmrz = 0.
        for iat in atoms:
            cmrx += st.atom[iat].x
            cmry += st.atom[iat].y
            cmrz += st.atom[iat].z
        n = len(atoms)
        cmrx /= float(n)
        cmry /= float(n)
        cmrz /= float(n)

        try:
            st = structure.StructureReader(self.ligidxfile).next()
        except:
            print "Cannot open Maestro structure file %s" % self.ligidxfile
            sys.exit(-1)
        nligatoms = len(st.atom)
        atoms = structureutil.evaluate_asl(st, ligand_asl)
        lig_atoms = [];
        for iat in atoms:
            lig_atoms.append(st.atom[iat].property['i_i_internal_atom_index']-nrecatoms)

        #computes ligand CM
        cmlx = cmly = cmlz = 0.
        for iat in atoms:
            cmlx += st.atom[iat].x
            cmly += st.atom[iat].y
            cmlz += st.atom[iat].z
        n = len(atoms)
        cmlx /= float(n)
        cmly /= float(n)
        cmlz /= float(n)

        #computes and reports the distance btw CM's
        d = sqrt((cmlx - cmrx)*(cmlx - cmrx) + (cmly - cmry)*(cmly - cmry) + (cmlz - cmrz)*(cmlz - cmrz))
        print "CM-CM distance = %f" % d

        self.restraint_file = self.jobname + '_cmrestraint.dat'
        f = open(self.restraint_file,"w")
        f.write("Receptor\n")
        f.write("%d\n" % 1)
        f.write("%d\n" % len(rec_atoms) )
        for i in rec_atoms:
            f.write("%d\n" % i)
        f.write("Ligand\n")
        f.write("%d\n" % 2)
        f.write("%d\n" % len(lig_atoms) )
        for i in lig_atoms:
            f.write("%d\n" % i)
        f.close()

#
# writes receptor mae file with restraints
#
    def writeRecStructureFile(self):
        receptor_file =  self.keywords.get('RECEPTOR_FILE')
        if not receptor_file:
            msg = "bedam_prep: No receptor file specified in the input file"
            self.exit(msg)
        if not os.path.exists(receptor_file):
            msg = 'File does not exist: %s' % receptor_file
            self.exit(msg)
        receptor_asl =  self.keywords.get('REST_RECEPTOR_ASL')

        try:
            st = structure.StructureReader(receptor_file).next()
        except:
            print "Cannot open Maestro structure file %s" % receptor_file
            sys.exit(-1)

        atoms = []
        if receptor_asl:
            atoms = structureutil.evaluate_asl(st, receptor_asl)
        
        # all atoms free
        for atom in st.atom:
            atom.property['i_i_constraint'] = 0
        # buffer the ones in ASL
        for iat in atoms:
            st.atom[iat].property['i_i_constraint'] = 2

        receptor_file_restr =   self.jobname + '_rcpt_restr' + '.maegz'
        st.write(receptor_file_restr)
        self.receptor_file_restr = receptor_file_restr

#
# writes the Impact input file for minimization/thermalization
#
    def  writeThermInputFile(self):
        if self.receptor_file_restr is None:
            msg = "writeThermInputFile: Internal error: receptor structure file not found"
            self.exit(msg)
        if not os.path.exists(self.receptor_file_restr):
            msg = 'File does not exist: %s' % self.receptor_file_restr
            self.exit(msg)
        ligand_file =  self.keywords.get('LIGAND_FILE')
        if not ligand_file:
            msg = "bedam_prep: No ligand file specified in the input file"
            self.exit(msg)
        if not os.path.exists(ligand_file):
            msg = 'File does not exist: %s' % ligand_file
            self.exit(msg)
        kfcm =  self.keywords.get('REST_LIGAND_CMKF')
        if not kfcm:
            kfcm = '3.0'
        d0cm =  self.keywords.get('REST_LIGAND_CMDIST0')
        if not d0cm:
            msg = "bedam_prep: No receptor-ligand reference distance specified"
            self.exit(msg)
        tolcm =  self.keywords.get('REST_LIGAND_CMTOL')
        if not tolcm:
            msg = "bedam_prep: No receptor-ligand distance tolerance specified"
            self.exit(msg)
        if self.restraint_file is None:
            msg = "writeThermInputFile: Internal error: restraint file not found"
            self.exit(msg)
        if not os.path.exists(self.restraint_file):
            msg = 'File does not exist: %s' % self.restraint_file
            self.exit(msg)
        temperature =  self.keywords.get('TEMPERATURE')
        if not temperature:
            temperature = '300.0'
        impact_input_file =   self.jobname + '_mintherm' + '.inp'
        impact_output_file =  self.jobname + '_mintherm' + '.out'
        impact_jobtitle =     self.jobname + '_mintherm'
        out_restart_file =    self.jobname + '_mintherm' + '.rst'
        out_structure_file =  self.jobname + '_mintherm' + '.maegz'
        self.mintherm_out_restart_file = out_restart_file 

        input = self.input_mintherm % (impact_output_file, impact_jobtitle, self.receptor_file_restr , ligand_file, kfcm, d0cm, tolcm, self.restraint_file, temperature, out_restart_file, out_structure_file)
        f = open(impact_input_file, "w")
        f.write(input)
        f.close


    def mkAuxParams(self):
        lambda_string = self.keywords.get('LAMBDAS')
        if not lambda_string:
            msg = "bedam_prep: No lambda schedule specified"
            self.exit(msg)
        lambdas = lambda_string.split(',')
        nlambdas = len(lambdas)
        qb_fcs_string = self.keywords.get('QB_FORCECONSTANT')
        if not qb_fcs_string:
            msg = "bedam_prep: No quadbias force constants specified"
            self.exit(msg)
        self.qb_fcs = qb_fcs_string.split(',')
        if not (len(self.qb_fcs) == nlambdas):
            msg = "bedam_prep: Number of quadbias force constants does not match the number of lambdas."
            self.exit(msg)
        qb_bepos_string = self.keywords.get('QB_BEPOS')
        if not qb_bepos_string:
            msg = "bedam_prep: No quadbias binding energy positions specified"
            self.exit(msg)
        self.qb_bepos = qb_bepos_string.split(',')
        if not (len(self.qb_bepos) == nlambdas):
            msg = "bedam_prep: Number of quadbias binding energy positions does not match the number of lambdas."
            self.exit(msg)


#
# writes the Impact input file for remd production
#
    def  writeRemdInputFile(self):
        if self.receptor_file_restr is None:
            msg = "writeRemdInputFile: Internal error: receptor structure file not found"
            self.exit(msg)
        if not os.path.exists(self.receptor_file_restr):
            msg = 'File does not exist: %s' % self.receptor_file_restr
            self.exit(msg)
        ligand_file =  self.keywords.get('LIGAND_FILE')
        if not ligand_file:
            msg = "bedam_prep: No ligand file specified in the input file"
            self.exit(msg)
        if not os.path.exists(ligand_file):
            msg = 'File does not exist: %s' % ligand_file
            self.exit(msg)
        umax = self.keywords.get('UMAX')
        if not umax:
            umax = '1000.0'
        kfcm =  self.keywords.get('REST_LIGAND_CMKF')
        if not kfcm:
            kfcm = '3.0'
        d0cm =  self.keywords.get('REST_LIGAND_CMDIST0')
        if not d0cm:
            msg = "bedam_prep: No receptor-ligand reference distance specified"
            self.exit(msg)
        tolcm =  self.keywords.get('REST_LIGAND_CMTOL')
        if not tolcm:
            msg = "bedam_prep: No receptor-ligand distance tolerance specified"
            self.exit(msg)
        if self.restraint_file is None:
            msg = "writeRemdInputFile: Internal error: restraint file not found"
            self.exit(msg)
        if not os.path.exists(self.restraint_file):
            msg = 'File does not exist: %s' % self.restraint_file
            self.exit(msg)
        temperature =  self.keywords.get('TEMPERATURE')
        if not temperature:
            msg = "bedam_prep: no temperature speciefied"
            self.exit(msg)
        

        lambda_string = self.keywords.get('LAMBDAS')
        if not lambda_string:
            msg = "bedam_prep: No lambda schedule specified"
            self.exit(msg)
        lambdas = lambda_string.split(',')
        nlambdas = len(lambdas)
        # constructs list of lambdas in input file
        lambda_list = ''
        for l in lambdas:
            lambda_list = lambda_list + "-\n      lambda %s " % l

        rest_kf = self.keywords.get('REST_RECEPTOR_KF')
        if not rest_kf:
            rest_kf = '0.6'

        if self.mintherm_out_restart_file is None:
            msg = "writeRemdInputFile: Internal error: min/therm restart file not found"
            self.exit(msg)
        basename = os.path.splitext( os.path.basename(self.mintherm_out_restart_file) )[0]
        suffix = os.path.splitext( os.path.basename(self.mintherm_out_restart_file) )[1]
        for i in range(0,nlambdas):
            os.system("ln -s %s %s_%d%s" % (self.mintherm_out_restart_file, basename, i, suffix))

        
        nmd_eq = self.keywords.get('EQUILIBRATION_STEPS')
        if not nmd_eq:
            msg = "bedam_prep: Number of equilibration steps not specified"
            self.exit(msg)
        nmd_prod = self.keywords.get('PRODUCTION_STEPS')
        if not nmd_prod:
            msg = "bedam_prep: Number of production steps not specified"
            self.exit(msg)
        nprnt = self.keywords.get('PRNT_FREQUENCY')
        if not nprnt:
            msg = "Number of printing frequency not specified"
            self.exit(msg)
        ntrj = self.keywords.get('TRJ_FREQUENCY')
        if not ntrj:
            msg = "Number of trajectory writing frequency not specified"
            self.exit(msg)

        rcpt_rsvr = self.keywords.get('RCPT_RESERVOIR_TRJFILE')
        lig_rsvr  = self.keywords.get('LIG_RESERVOIR_TRJFILE')
        nrsvr     = self.keywords.get('RESERVOIR_SIZE')
        if rcpt_rsvr and lig_rsvr and nrsvr:
            if not os.path.exists(rcpt_rsvr):
                msg = 'File does not exist: %s' % rcpt_rsvr
                self.exit(msg)
            if not os.path.exists(rcpt_rsvr):
                msg = 'File does not exist: %s' % lig_rsvr
                self.exit(msg)
            rsvr_command = self.input_reservoir % (nrsvr, nrsvr, rcpt_rsvr, lig_rsvr)
        else:
            rsvr_command = ""

        restart_file = self.keywords.get('RESTART_FILE')

        #Bias potentials settings
        bias_flag = ""
        bias_parameters = "! "
        qb_k = self.keywords.get('QUADBIAS')
        if qb_k:
            if qb_k == "yes":
                bias_flag = "quadbias"
                self.mkAuxParams() #creates self.qb_fcs and self.qb_bepos
                # constructs list of auxiliary bedam parameters in input file
                n_aux_parameters = 2 # force constant + binding energy position
                bias_parameters = "    weight auxbind naux %d " % n_aux_parameters
                for i in range(0,nlambdas):
                    bias_parameters = bias_parameters + "-\n      auxp1 %s auxp2 %s " % (self.qb_fcs[i],self.qb_bepos[i])

        impact_input_file =   self.jobname + '_remd' + '.inp'
        impact_output_file =  self.jobname + '_remd' + '.out'
        impact_jobtitle =     self.jobname + '_remd'
        out_trajectory_file = self.jobname + '_remd' + '.trj'
        out_restart_file =    self.jobname + '_remd' + '.rst'
        out_structure_file =  self.jobname + '_remd' + '.maegz'


        if not restart_file:
            input = self.input_remd % (impact_output_file, impact_jobtitle, 
                                   self.receptor_file_restr , ligand_file,
                                   umax, bias_flag,
                                   rest_kf, nlambdas, lambda_list,
                                   bias_parameters,
                                   kfcm, d0cm, tolcm, self.restraint_file, 
                                   self.mintherm_out_restart_file,
                                   nmd_eq, nlambdas, nmd_eq, temperature, temperature, nprnt,
                                   nmd_prod, nlambdas, rsvr_command,
                                   temperature, temperature, nprnt,
                                   ntrj, out_trajectory_file, out_restart_file, out_structure_file)
        else:
            input = self.input_remd_restart % (impact_output_file, impact_jobtitle, 
                                   self.receptor_file_restr , ligand_file,
                                   umax, bias_flag,
                                   rest_kf, nlambdas, lambda_list,
                                   bias_parameters,
                                   kfcm, d0cm, tolcm, self.restraint_file, 
                                   nmd_prod, nlambdas, rsvr_command,
                                   temperature, temperature, nprnt,
                                   restart_file,
                                   ntrj, out_trajectory_file, 
                                   out_restart_file, out_structure_file)

        f = open(impact_input_file, "w")
        f.write(input)
        f.close



#
# writes the Impact input file to construct .maegz trajectories by lambda
#
    def  writeReadTrajInputFile(self):
        if self.receptor_file_restr is None:
            msg = "writeReadTrajInputFile: Internal error: receptor structure file not found"
            self.exit(msg)
        if not os.path.exists(self.receptor_file_restr):
            msg = 'File does not exist: %s' % self.receptor_file_restr
            self.exit(msg)
        ligand_file =  self.keywords.get('LIGAND_FILE')
        if not ligand_file:
            msg = "bedam_prep: No ligand file specified in the input file"
            self.exit(msg)
        if not os.path.exists(ligand_file):
            msg = 'File does not exist: %s' % ligand_file
            self.exit(msg)
        umax = self.keywords.get('UMAX')
        if not umax:
            umax = '1000.0'
        lambda_string = self.keywords.get('LAMBDAS')
        if not lambda_string:
            msg = "bedam_prep: No lambda schedule specified"
            self.exit(msg)
        lambdas = lambda_string.split(',')
        nlambdas = len(lambdas)
        # constructs list of lambdas in input file
        lambda_list = ''
        for l in lambdas:
            lambda_list = lambda_list + "-\nlambda %s " % l

        rest_kf = self.keywords.get('REST_RECEPTOR_KF')
        if not rest_kf:
            rest_kf = '0.6'

        impact_input_file =   self.jobname + '_readtraj' + '.inp'
        impact_output_file =  self.jobname + '_readtraj' + '.out'
        impact_jobtitle =     self.jobname + '_readtraj'

        # creates short-named versions of trj and maegz files
        for i in range(nlambdas):
            long_trjfile = self.jobname + "_remd_%d.trj" % i
            short_trjfile = "q%d.trj" % i
            os.system("ln -s -f %s %s" % (long_trjfile, short_trjfile))

        
        i = 0
        for l in lambdas:
            long_maefile = self.jobname + "_remd_trj_%s.maegz" % l
            short_maefile = "q%d.maegz" % i
            os.system("ln -s -f %s %s" % (short_maefile, long_maefile))
            i += 1

        input = self.input_trj % (impact_output_file, impact_jobtitle, 
                                   self.receptor_file_restr , ligand_file, 
                                  umax,
                                   rest_kf, nlambdas, lambda_list, 
                                   nlambdas)
        f = open(impact_input_file, "w")
        f.write(input)
        f.close

        #writes readtraj for surface rescoring
        impact_input_file =   self.jobname + '_readtraj_rescore' + '.inp'
        impact_output_file =  self.jobname + '_readtraj_rescore' + '.out'
        impact_jobtitle =     self.jobname + '_readtraj_rescore'
        input = self.input_trj_rescore % (impact_output_file, impact_jobtitle, 
                                   self.receptor_file_restr , ligand_file, 
                                  umax,
                                   rest_kf, nlambdas, lambda_list, 
                                   nlambdas,nlambdas-1)
        f = open(impact_input_file, "w")
        f.write(input)
        f.close


#
# Reads all of the simulation data values temperature, energies, etc.
# at each time step and puts into a big table
#
    def getImpactData(self,file):
        if not os.path.exists(file):
            msg = 'File does not exist: %s' % file
            self.exit(msg)
        step_line = re.compile("^ Step number:")
        number_line = re.compile("(\s+-*\d\.\d+E[\+-]\d+\s*)+")
        temperature_line = re.compile("^\s*input target temperature\s+(\d*\.*\d*)")
        have_trgtemperature = 0
        nsamples = 0
        data = []
        f = open(file ,"r")
        line = f.readline()
        while line:
            # fast forward until we get to the line: 
            # "Step number: ... ", grab target temperature along the way
            # if it's in the input file
            while line and not re.match(step_line, line):
                if re.match(temperature_line, line):
                    words = line.split()
                    temperature = words[3]
                    have_trgtemperature = 1
                line = f.readline()
            # read the step number
            if re.match(step_line, line):
                words = line.split()
                step = words[2]
                datablock = [int(step)]
                if have_trgtemperature == 1:
                    datablock.append(float(temperature))
                #now read up to 3 lines of numbers
                ln = 0
                while ln < 3:
                    line = f.readline()
                    if not line:
                        msg = "Unexpected end of file"
                        self.exit(msg)
                    if re.match(number_line, line):
                        for word in line.split():
                            datablock.append(float(word))
                        ln += 1
                data.append(datablock)
            line = f.readline()
        f.close()
        return data
        """
            #read up to 3 lines of numbers
            ln = 0
            datablock = []
            while ln < 3:
                line = f.readline()
                if re.match(number_line,line):
                    print line
                    ln += 1
            line = f.readline()
        """
                
        """
            # fast forward until we get to the line: 
            # "Step number: ... "
            for line in f:
                # fast forward until we get to the line: 
                # "Step number: ... "
                if not in_block:
                    if not (line.find("Step number:") == -1):
                        in_block = 1
                        words = line.split()

                line = f.readline()
            while(line.find('Step number:') == -1):
        """ 


#
# Extract binding energies from data table
#
    def getBindingEnergies(self):
        """
        extract binding energies from .out files
        """
        lambda_string = self.keywords.get('LAMBDAS')
        if not lambda_string:
            msg = "getBindingEnergies: No lambda schedule specified"
            self.exit(msg)
        lambdas = lambda_string.split(',')
        nlambdas = len(lambdas)
        self.lambdas = []
        for i in range(0,nlambdas):
            self.lambdas.append(float(lambdas[i]))
            
        nmd_eq = self.keywords.get('EQUILIBRATION_STEPS')
        if not nmd_eq:
            msg = "bedam_prep: Number of equilibration steps not specified"
            self.exit(msg)
        nmd_prod = self.keywords.get('PRODUCTION_STEPS')
        if not nmd_prod:
            msg = "bedam_prep: Number of production steps not specified"
            self.exit(msg)
        nprnt = self.keywords.get('PRNT_FREQUENCY')
        if not nprnt:
            msg = "Number of printing frequency not specified"
            self.exit(msg)
        nreject = int(nmd_eq)/int(nprnt)
        
        self.binding_energies_bysim = []
        self.lambda_bysim = []
        self.binding_energies_bylambda = {}
        for i in range(0,nlambdas):
            impact_output_file =  self.jobname + '_remd' '_%d' % i + '.out'
            datai = self.getImpactData(impact_output_file)
            be = []
            lv = []
            ndata = len(datai)
            nf = len(datai[0])
            for i in range(nreject,ndata):
                lmb = datai[i][nf-2]
                u = datai[i][nf-1]
                be.append(u)
                lv.append(lmb)
                if not lmb in self.binding_energies_bylambda.keys():
                    self.binding_energies_bylambda[lmb] = []
                self.binding_energies_bylambda[lmb].append(u)
            self.binding_energies_bysim.append(be)
            self.lambda_bysim.append(lv)



# Return the indexes of the output files to be analyzed
    def collect_production_cycles2D(self,replica,cycles_eq,basename):
        replica_dir = "r%d" % replica
        out_files = glob.glob(replica_dir + '/' + basename + '_*.out')
	to_cycle = re.compile(replica_dir + '/' + basename + r"_(\d+).out")
	cycles = []
        for f in out_files:
		c = re.match(to_cycle, f).group(1)
		cycles.append(int(c))
	cycles.sort()
	return cycles[cycles_eq:]

#
# Extract binding energies and total energies from data table (for 2D T/lambda RE)
#
    def getTotalAndBindingEnergies(self):
        """
        extract binding energies and total energies from .out files
        """
        lambda_string = self.keywords.get('LAMBDAS')
        if not lambda_string:
            msg = "getTotalAndBindingEnergies: No lambda schedule specified"
            self.exit(msg)
        lambdas = lambda_string.split(',')
        if self.keywords.get('TEMPERATURES') is None:
            self._exit("getTotalAndBindingEnergies: No temperature schedule specified")
        temperatures = self.keywords.get('TEMPERATURES').split(',')
        #build parameters for the lambda/temperatures combined states
        self.nreplicas = self._buildBEDAMStates(lambdas,temperatures)

        print self.nreplicas
        
        nlambdas = len(lambdas)        
        self.lambdas = []
        for i in range(0,nlambdas):
            self.lambdas.append(float(lambdas[i]))
        ntemperatures = len(temperatures)
        self.temperatures = []
        for i in range(0,ntemperatures):
            self.temperatures.append(float(temperatures[i]))

        basename = self.keywords.get('ENGINE_INPUT_BASENAME')
        if not basename:
            msg = "No job basename (ENGINE_INPUT_BASENAME) specified"
            self.exit(msg)

        cycles_eq = self.keywords.get('EQUILIBRATION_CYCLES')
        if not cycles_eq:
            msg = "bedam_analyze (2D): Number of equilibration cycles not specified"
            self.exit(msg)
        cycles_eq = int(cycles_eq)

        # bias info contains:
        # (temperature, total_energy, lambda, binding_energy)
        # ordered by simulation time (within a simulation)
        self.bias_info_bysim = [[] for i in range (0,self.nreplicas)]
        self.nsamples_sim = [0 for i in range(0,self.nreplicas)]
        self.nsamples_state = {}
        self.bias_info_bystate = {}
        for repl in range(0,self.nreplicas):
            self.nsamples_sim[repl] = 0
            cycles = self.collect_production_cycles2D(repl,cycles_eq, basename)
            for cy in cycles:
                impact_output_file =  'r%d' % repl + '/' + basename + '_%d' % cy + '.out'
                try:
                    datai = self.getImpactData(impact_output_file)
                except:
                    print "Warning: problems reading file %s" % impact_output_file
                    continue
                ndata = len(datai)
                if ndata == 0:
                    print "Warning: problems reading file %s" % impact_output_file
                    continue
                nf = len(datai[0])
                for i in range(0,ndata):
#                for i in range(ndata-1,ndata):
                    lmb = datai[i][nf-2]
                    u = datai[i][nf-1]
                    temperature = datai[i][1]
                    total_energy = datai[i][3]
                    info = {
                            "temperature": temperature,
                            "total_energy": total_energy,
                            "lambda": lmb,
                            "binding_energy": u
                    }
                    self.bias_info_bysim[repl].append(info)
                    self.nsamples_sim[repl] += 1;
                    k = (temperature, lmb)
                    if not k in self.nsamples_state.keys():
                        self.nsamples_state[k] = 0 
                        self.bias_info_bystate[k] = []
                    self.nsamples_state[k] += 1
                    self.bias_info_bystate[k].append(info)

    def _buildBEDAMStates(self,lambdas,temperatures):
        self.stateparams = []
        for lambd in lambdas:
            for tempt in temperatures:
                st = {}
                st['lambda'] = lambd
                st['temperature'] = tempt
                self.stateparams.append(st)
        return len(self.stateparams)

#
# Creates a binding energy grid
#
    def mkdefault_grid(self):
        
        grid_tail = [2,4,8,16,32,64,128,256,512,999.999,1001]
        dbe = 0.5
        # find minimum energy
        mk = []
        for k in range(len(self.binding_energies_bysim)):
            mk.append(min(array(self.binding_energies_bysim[k],dtype=float)))
        tm = min(mk) - dbe
        ugrid = arange(start=tm,stop=1.,step=dbe,dtype=float)
        ugrid = append(ugrid,array(grid_tail,dtype=float))
        return ugrid

#
# Saves histograms of binding energies in a/be_hist/h_<lambda>.dat
#
    def be_histograms(self):
        if self.binding_energies_bylambda is None:
            self.getBindingEnergies()

        grid_string = self.keywords.get('BE_GRID')
        if not grid_string:
            ugrid = self.mkdefault_grid()
        else:
            # numpy array
            ugrid = array(grid_string.split(','),dtype=float)

        os.system("mkdir -p a/be_hist")
        for lmb in self.binding_energies_bylambda.keys():
            # numpy histogram function
            h,bins = histogram(array(self.binding_energies_bylambda[lmb]),bins=ugrid)
            file = "a/be_hist/h_%f.dat" % lmb;
            f = open(file ,"w")
            for i in range(len(bins)-1):
                f.write("%14.4e %d\n" % (bins[i], h[i]))
            f.close

#
# Saves be/lambda trajectories for each simulation in a/lbe_trj/lbe_<sim>.dat
#
    def lbe_trj(self):
        if self.binding_energies_bysim is None:
            self.getBindingEnergies()

        os.system("mkdir -p a/lbe_trj")
        for k in range(len(self.binding_energies_bysim)):
            file = "a/lbe_trj/lbe_%d.dat" % k;
            f = open(file ,"w")
            for i in range(len(self.binding_energies_bysim[k])):
                f.write("%f %f\n" % (self.lambda_bysim[k][i],self.binding_energies_bysim[k][i]))
            f.close 
        

#
# Saves samples (not histograms!) of binding energies/total energies in a/be_hist/h_<i>.dat
#
    def be_histograms_lt(self):
        if self.bias_info_bystate is None:
            self.getTotalAndBindingEnergies()

        os.system("mkdir -p a/be_hist")
        for l in sorted(self.bias_info_bystate.keys(), key=itemgetter(1,0)):
            (tempt,lmb) = l
            file = "a/be_hist/h_%f_%f.dat" % (lmb,tempt);
            f = open(file ,"w")
            for entry in self.bias_info_bystate[l]:
                u = entry['binding_energy']
                e = entry['total_energy']
                f.write("%16.6e %16.6e\n" % (u,e))
            f.close

#
# Saves be/total energy/lambda/temperature trajectories for each replica in a/lbe_trj/lbe_<sim>.dat
#
    def lbe_trj_lt(self):
        if self.bias_info_bysim is None:
            self.getTotalAndBindingEnergies()

        os.system("mkdir -p a/lbe_trj")
        for k in range(0,self.nreplicas):
            file = "a/lbe_trj/lbe_%d.dat" % k;
            f = open(file ,"w")
            for info in self.bias_info_bysim[k]:             
                f.write("%f %f %f %f\n" % (info["temperature"],info["total_energy"],info["lambda"],info["binding_energy"]))
            f.close 







#
# Runs MBAR using binding energies
#
    def runMBAR(self):
        """
        runs MBAR using lambdas and the collected binding energies
        """
        if self.lambdas is None:
            msg = "runMBAR: Internal error: lambdas are not defined"
            self.exit(msg)
        if self.binding_energies_bysim is None:
            msg = "runMBAR: Internal error: binding energies are not defined"
            self.exit(msg)

        K = len(self.lambdas)
        n = len(self.binding_energies_bysim[0])
        N_k = zeros([K], int32)
        for i in range(K):
            N_k[i] = n

        temperature =  float(self.keywords.get('TEMPERATURE'))        
        if not temperature:
            msg = "bedam_prep: no temperature specified"
            self.exit(msg)
        beta = 1.0 / (0.001986209 * temperature)

        #compute Vsite term
        kfcm =  float(self.keywords.get('REST_LIGAND_CMKF'))
        if not kfcm:
            kfcm = '3.0'
        d0cm =  self.keywords.get('REST_LIGAND_CMDIST0')
        if not d0cm:
            msg = "bedam: No receptor-ligand reference distance specified"
            self.exit(msg)
        d0cm = float(d0cm)
        tolcm =  float(self.keywords.get('REST_LIGAND_CMTOL'))
        if not tolcm:
            msg = "bedam: No receptor-ligand distance tolerance specified"
            self.exit(msg)
        r1 =  d0cm - tolcm
        if r1 < 0.:
            r1 = 0.0
        r2 = d0cm + tolcm
        pi = acos(-1.)
        vsite = (4./3.)*pi*(r2*r2*r2 - r1*r1*r1)
        c0 = 1./1660.
        dgvsite = -(1./beta)*log(c0*vsite)
        print "dgvsite = %f" % dgvsite

        #run MBAR
        u_klt = zeros([K,K,n], float64)
        u_kt = zeros([K,n], float64)
        for k in range(K):
            for t in range(N_k[k]):
                u_kt[k,t] = self.binding_energies_bysim[k][t]
        for k in range(K):
            for t in range(N_k[k]):
                for l in range(K):
                    # Reduced energy of snapshot t from simulation k 
                    # in state l
                    u_klt[k,l,t] = beta * self.lambdas[l] * u_kt[k,t]

        qb_k = self.keywords.get('QUADBIAS')
        if qb_k:
            if qb_k == "yes":
                self.mkAuxParams()
                for k in range(K):
                    for t in range(N_k[k]):
                        for l in range(K):
                            # add quadratic bias
                            du = u_kt[k,t]-float(self.qb_bepos[l])
                            u_klt[k,l,t] += beta * 0.5 * float(self.qb_fcs[l]) * du*du

        print "Running MBAR..."
        mbar = MBAR.MBAR(u_klt, N_k, relative_tolerance = 1.0e-6, verbose = True )
        (Deltaf_ij, dDeltaf_ij) = mbar.getFreeEnergyDifferences()

        #dimensionless free energies
        fk = zeros([K],float64)
        for i in range(len(Deltaf_ij)):
            fk[i] = Deltaf_ij[0][i]
        print fk

        #weight of each sample
        weights_bysim = zeros([K,n], float64)
        wt = zeros([K], float64)
        for k in range(K):
            for t in range(N_k[k]):
                # Reduced energy of snapshot t from simulation k 
                # in umbrella potential 
                wt = N_k[:]*exp(fk[:] - u_klt[k,:,t])
                weights_bysim[k][t] = 1./sum(wt)
        weights_bysim = weights_bysim/sum(weights_bysim)

        #computes p0
        grid_string = self.keywords.get('BE_GRID')
        if not grid_string:
            ugrid = self.mkdefault_grid()
        else:
            # numpy array
            ugrid = array(grid_string.split(','),dtype=float64)

        h,bins = histogram(u_kt,bins=ugrid,weights=weights_bysim)
        file = "a/be_hist/p0.dat";
        f = open(file ,"w")
        for i in range(len(bins)-1):
            f.write("%16.6e %16.6e %16.6e\n" % (bins[i], h[i], h[i]/(bins[i+1]-bins[i])))
        f.close

        #convert to kcal/mol
        for i in range(len(Deltaf_ij)):
            for j in range(len(Deltaf_ij[i])):
                Deltaf_ij[i][j] = Deltaf_ij[i][j]/beta
                dDeltaf_ij[i][j] = dDeltaf_ij[i][j]/beta
        #add Vsite term
        for i in range(len(Deltaf_ij)):
            for j in range(len(Deltaf_ij[i])):
                Deltaf_ij[i][j] += dgvsite

        return (Deltaf_ij, dDeltaf_ij)


#
# Runs MBAR using binding energies
#
    def runMBAR_BEDAMTempt(self):
        """
        runs MBAR using lambdas and the collected binding energies
        """
        if self.lambdas is None:
            msg = "runMBAR: Internal error: lambdas are not defined"
            self.exit(msg)
        if self.temperatures is None:
            msg = "runMBAR: Internal error: lambdas are not defined"
            self.exit(msg)
        if self.bias_info_bysim is None:
            msg = "runMBAR: Internal error: bias info is not available"
            self.exit(msg)

        K = self.nreplicas
        N_k = zeros([K], int32)
        
        # states are defined by (temperature, lambda) tuples. This orders states by lambda first
        # from 0 to 1 and then by temperature
        states = sorted(self.nsamples_state.keys(), key=itemgetter(1,0))

        i = 0
        nmax = 0
        for l in states:
            N_k[i] = self.nsamples_state[l]
            if N_k[i] > nmax:
                nmax = N_k[i]
            i += 1
        if i != K:
            msg = "runMBAR_BEDAMTempt: internal error, inconsistency in number of states."
            self.exit(msg)

        lambda_ref = 0.0
        temperature_ref =  float(self.keywords.get('TEMPERATURE'))        
        if not temperature_ref:
            msg = "bedam_prep: no temperature specified"
            self.exit(msg)
        beta_ref = 1.0 / (0.001986209 * temperature_ref)

        beta = {}
        for l in states:
            (tempt,lmb) = l
            beta[l] = 1.0 / (0.001986209 * tempt)

        #run MBAR
        u_klt = zeros([K,K,nmax], float64)
        for k in range(K):
            for t in range(self.nsamples_sim[k]):
                u = self.bias_info_bysim[k][t]['binding_energy']
                e = self.bias_info_bysim[k][t]['total_energy']
                li = self.bias_info_bysim[k][t]['lambda']
                # this is H0(x) = H(x) - lambda u(x), here li is the lambda from the
                # simulation when the sample was collected. (This is because the total energy
                # includes lambda*u )
                e0 = e - li*u
                i = 0
                for l in states:
                    (tempt,lmb) = l
                    b = beta[l]
                    # beta_i * lambda_i - beta_ref * lambda_ref
                    bl = b * lmb - beta_ref * lambda_ref
                    # Reduced energy of snapshot t from simulation k 
                    # in state l
                    u_klt[k,i,t] = (b - beta_ref) * e0 + bl * u
                    i += 1

        print "Running MBAR..."
        mbar = MBAR.MBAR(u_klt, N_k, maximumIterations = 1000, relative_tolerance = 1.0e-6, verbose = True )
        (Deltaf_ij, dDeltaf_ij) = mbar.getFreeEnergyDifferences(uncertainty_method = 'none')

        #dimensionless free energies
        fk = zeros([K],float64)
        for i in range(len(Deltaf_ij)):
            fk[i] = Deltaf_ij[0][i]
        print fk

        #weight of each sample
#        weights_bysim = zeros([K,n], float64)
#        wt = zeros([K], float64)
#        for k in range(K):
#            for t in range(N_k[k]):
                # Reduced energy of snapshot t from simulation k 
                # in umbrella potential 
#                wt = N_k[:]*exp(fk[:] - u_klt[k,:,t])
#                weights_bysim[k][t] = 1./sum(wt)
#        weights_bysim = weights_bysim/sum(weights_bysim)

        #computes p0
#        grid_string = self.keywords.get('BE_GRID')
#        if not grid_string:
#            ugrid = self.mkdefault_grid()
#        else:
            # numpy array
#            ugrid = array(grid_string.split(','),dtype=float64)

#       h,bins = histogram(u_kt,bins=ugrid,weights=weights_bysim)
#       file = "a/be_hist/p0.dat";
#       f = open(file ,"w")
#       for i in range(len(bins)-1):
#           f.write("%16.6e %16.6e %16.6e\n" % (bins[i], h[i], h[i]/(bins[i+1]-bins[i])))
#       f.close

        #convert to kcal/mol
#        for i in range(len(Deltaf_ij)):
#            for j in range(len(Deltaf_ij[i])):
#                Deltaf_ij[i][j] = Deltaf_ij[i][j]/beta_ref
#                dDeltaf_ij[i][j] = dDeltaf_ij[i][j]/beta_ref
        #add Vsite term
#        for i in range(len(Deltaf_ij)):
#            for j in range(len(Deltaf_ij[i])):
#                Deltaf_ij[i][j] += dgvsite

        return (Deltaf_ij, dDeltaf_ij)

    
