# bedam python class
__doc__="""
$Revision: 0.1 $

A class to prepare BEDAM AsyncRE jobs

"""
# Contributors: Emilio Gallicchio

import os, sys, time, re, glob
from schrodinger.utils import cmdline
import schrodinger.utils.log
import shutil
import signal
import glob
import time

import sqlite3 as lite

from math import *

from bedam_prep_ac import bedam_prep_ac

class bedam_job_asyncre(bedam_prep_ac):
    """
    Class to set up BEDAM calculations
    """
    def __init__(self, command_file, options):
        bedam_prep_ac.__init__(self, command_file, options)

#
# Impact input file templates for academic
#
    def setupTemplatesASyncRE(self):
        """ Setup templates for input files for academic impact"""
        self.input_cms =  """
task {
  task = "desmond:auto"
}

build_geometry {
  box = {
     shape = "orthorhombic"
     size = [10.0 10.0 10.0 ]
     size_type = "absolute"
  }
  neutralize_system = false
  rezero_system = false
  solvate_system = false
}

assign_forcefield {
}

"""  

        self.input_idx =  """
write file -
"{out_file}" -
      title -
"{title}" *

CREATE
  build primary name species1 type auto read sqldb file -
"{dms_rcpt_in}"
  build primary name species2 type auto read sqldb file -
"{dms_lig_in}"
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
  input cntl mxcyc 0  rmscut 0.05 deltae 1.0e-05
  conjugate dx0 0.05 dxm 1.0
  run
  write sql name species1 file "{dms_rcpt_out}"
  write sql name species2 file "{dms_lig_out}"
QUIT

END
"""
        self.input_mintherm = """
write file -
"{out_file}" -
      title -
"{title}" *

CREATE
  build primary name species1 type auto read sqldb file -
"{dms_rcpt_in}"
  build primary name species2 type auto read sqldb file -
"{dms_lig_in}"
QUIT

SETMODEL
  setpotential
    mmechanics consolv agbnp2
  quit
  read parm file -
"paramstd.dat" -
  noprint
  energy rest domain cmdist kdist {cmkf} dist0 {cmdist0} toldist {cmtol} -
      read file "{cmrestraints_file}"
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
put {temperature} into 'tempt'
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
  write restart coordinates formatted file "{rst_file_out}"  
  write sql name species1 file "{dms_rcpt_out}"
  write sql name species2 file "{dms_lig_out}"
QUIT

END
"""
        self.input_remd = """
write file -
"{job_name}_@n@.out" -
      title -
"{job_name}" *

CREATE
  build primary name species1 type auto read sqldb file -
"{dms_rcpt_in}"
  build primary name species2 type auto read sqldb file -
"{dms_lig_in}"
QUIT

SETMODEL
  setpotential
    mmechanics nb12softcore umax {umax0} consolv agbnp2
    weight constraints buffer {rest_kf}
    weight bind rxid 0 nrep 1 lambda -
@lambda@
  quit
  read parm file -
"paramstd.dat" -
  noprint
  energy rest domain cmdist kdist {cmkf} dist0 {cmdist0} toldist {cmtol} -
      read file "{cmrestraints_file}"
  energy parm dielectric 1 nodist -
   listupdate 10 -
    cutoff 12 hmass 5
  energy rescutoff byatom all
  zonecons auto
  energy constraints bonds hydrogens
QUIT

if @n@ eq 1
DYNAMICS
  read restart coordinates formatted file "{job_name}_@nm1@.rst"
  input target temperature @temperature@
  input cntl initialize temperature at @temperature@
QUIT
endif

if @n@ gt 1
DYNAMICS
  read restart coordinates and velocities formatted file "{job_name}_@nm1@.rst"
QUIT
endif

!equilibration at new state
DYNAMICS
  input cntl nstep {nmd_eq} delt 0.0005
  input cntl constant temperature langevin relax 0.25
  input target temperature @temperature@
  input cntl nprnt {nmd_eq}
  input cntl tol 1.00000e-07
  input cntl stop rotations
  input cntl statistics off
  run rrespa fast 8
QUIT

DYNAMICS
  input cntl nstep {nmd} delt 0.001
  input cntl constant temperature langevin relax 1.0
  input target temperature @temperature@
  input cntl nprnt {nprnt}
  input cntl tol 1.00000e-07
  input cntl stop rotations
  input cntl statistics off
  run rrespa fast 4
  write restart coordinates and velocities formatted file "{job_name}_@n@.rst"
  write sql file "{job_name}_rcpt_@n@.dms" name species1
  write sql file "{job_name}_lig_@n@.dms" name species2
QUIT


END
"""

        self.input_slurm = """
#!/bin/bash
#SBATCH -p normal-mic       # Queue name
#SBATCH -N 6                # This is nodes, not cores (16 cores/node)
#SBATCH -n 6                # one process per node so we get one entry per node
#SBATCH -t 47:00:00         # Max time allotted for job
#SBATCH -A TG-MCB100145

echo "Number of nodes: $SLURM_NNODES"
echo "Nodelist: $SLURM_NODELIST"
echo "Number of tasks: $SLURM_NTASKS"
echo "Tasks per node: $SLURM_TASKS_PER_NODE"

scontrol show hostname $SLURM_NODELIST > .slurm_nodes
awk '{{ for(i=0;i<4;i++)print $1 ","i",4,Linux-x86_64,,/tmp"}}; {{ for(i=0;i<10;i++
)print $1"-mic0,"i",24,Linux-mic,,/tmp"}}; ' < .slurm_nodes > nodefile

python ~/src/async_re-0.3.2-alpha-multiarch/bedamtempt_async_re.py {job_name}_asyncre.cntl > LOG 2>&1
"""


        self.input_qsub = """
#!/bin/bash
#PBS -q production
#PBS -l select=64:ncpus=1
#PBS -N {job_name}
#PBS -l place=free
#PBS -l walltime=48:10:00
#PBS -V

cd $PBS_O_WORKDIR
sdir=/scratch/e.gallicchio

source ${{sdir}}/env/bin/activate

cp $PBS_NODEFILE .qsub_nodes
#1 core per replica
awk '{{ for(i=0;i<1;i++)print $1 ","i",1,Linux-x86_64,,/tmp"}}' < .qsub_nodes > nodefile

python ~/src/async_re-0.3.2-alpha-multiarch/bedamtempt_async_re.py {job_name}_asyncre.cntl > LOG 2>&1
"""

        self.input_runimpact_standard = """#!/bin/bash
export IMPACTHOME={acd_impact_home}
export IMP_ROOT=$IMPACTHOME
export IMPACT_EXEC=$IMP_ROOT/bin/Linux-x86_64
export LD_LIBRARY_PATH=$IMP_ROOT/lib/Linux-x86_64:$LD_LIBRARY_PATH
export OMP_NUM_THREADS={subjob_cores}
nice $IMPACT_EXEC/main1m $1
"""

        self.input_runimpact_multiarch = """#!/bin/bash
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
nice ./main1m $1
"""

        self.input_runimpact_interactive = """#!/bin/bash
export IMP_ROOT=$IMPACTHOME
export IMPACT_EXEC=$IMP_ROOT/bin/Linux-x86_64
export LD_LIBRARY_PATH=$IMP_ROOT/lib/Linux-x86_64:$LD_LIBRARY_PATH
$IMPACT_EXEC/main1m $1
"""


#
# Convert .mae files into .dms files with AGBNP2 parameters and internal atom indexes
#
    def getDesmondDMSFiles(self):
 
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
        com_shrod_source =  self.keywords.get('COMMERCIAL_SCHRODINGER_EVN')
        if not com_shrod_source:
            msg = "bedam_prep: No commerical shrodinger source file specified in the input file"
            self.exit(msg)

        print "Convert maegz files to cms files "
        desmond_builder_file = 'des_builder.msj'
        rcpt_cms_file = self.jobname + '_rcpt.cms'
        lig_cms_file = self.jobname + '_lig.cms'
        f = open(desmond_builder_file, 'w')
        input =  self.input_cms
        f.write(input)
        f.close()
        source_cmd = '. ' + com_shrod_source
        rcpt_cmd = '$SCHRODINGER/utilities/multisim' + ' -JOBNAME ' + self.jobname + ' -m ' + desmond_builder_file + ' ' + receptor_file + ' -o ' + rcpt_cms_file + ' -HOST localhost -maxjob 1 -WAIT'
        lig_cmd =  '$SCHRODINGER/utilities/multisim' + ' -JOBNAME ' + self.jobname + ' -m ' + desmond_builder_file + ' ' + ligand_file + ' -o ' + lig_cms_file + ' -HOST localhost -maxjob 1 -WAIT' 
        cms_cmd = source_cmd + ";" + rcpt_cmd + ";" + lig_cmd 
        os.system(cms_cmd)

        print "Convert cms files to dms files"
        rcpt_dms_file = self.jobname + '_rcpt.dms'
        lig_dms_file = self.jobname + '_lig.dms'
	rcpt_cmd = '$SCHRODINGER/run -FROM desmond mae2dms ' + rcpt_cms_file + ' ' + rcpt_dms_file
        lig_cmd = '$SCHRODINGER/run -FROM desmond mae2dms ' + lig_cms_file + ' ' + lig_dms_file
        dms_cmd = source_cmd + ";" + rcpt_cmd + ";" + lig_cmd
        os.system(dms_cmd)

        print "add agbnp parameters into dms files"
        agbnp_cmd =  "$SCHRODINGER/run add_agbnp2.py " + rcpt_dms_file
        agbnp_cmd = source_cmd + ";" + agbnp_cmd
        os.system(agbnp_cmd)
        agbnp_cmd =  "$SCHRODINGER/run add_agbnp2.py " + lig_dms_file
        agbnp_cmd = source_cmd + ";" + agbnp_cmd
        os.system(agbnp_cmd)

        print "add internal atom indexes into dms files"
        acd_impact_source =  self.keywords.get('ACADEMIC_IMPACT_EVN')
        if not acd_impact_source:
            msg = "bedam_prep: No academic IMPACT source file specified in the input file"
            self.exit(msg)
        impact_input_file =   self.jobname + '_idx' + '.inp'
        impact_output_file =  self.jobname + '_idx' + '.out'
        impact_jobtitle =     self.jobname + '_idx'
        out_receptor_file =   self.jobname + '_rcpt_idx' + '.dms'
        out_ligand_file =     self.jobname + '_lig_idx' + '.dms'
        f = open(impact_input_file, 'w')
        input =  self.input_idx.format(out_file=impact_output_file, title=impact_jobtitle, dms_rcpt_in=rcpt_dms_file, dms_lig_in=lig_dms_file, dms_rcpt_out=out_receptor_file, dms_lig_out=out_ligand_file )
        f.write(input)
        f.close()
        idx_log_file =  self.jobname + '_idx' + '.log'
        source_cmd = '. ' + acd_impact_source
        idx_cmd =  "$IMPACT_EXEC/main1m " + impact_input_file + " > " + idx_log_file + " 2>&1 "
        idx_cmd =  source_cmd + ";" + idx_cmd
        os.system(idx_cmd)
        if not os.path.exists(out_receptor_file) or not os.path.exists(out_ligand_file):
            print " failed"
            msg = "Impact job to generate dms files with idx failed"
            self.exit(msg)
        os.rename(out_receptor_file, rcpt_dms_file)
        os.rename(out_ligand_file, lig_dms_file)
        self.recidxfile = rcpt_dms_file 
        self.ligidxfile = lig_dms_file

    def writeCntlFile(self):
        input = ""

        job_transport = self.keywords.get('JOB_TRANSPORT')
        if job_transport is None:
            msg = "writeCntlFile: JOB_TRANSPORT is not specified"
            self.exit(msg)
        if not (job_transport == "SSH" or job_transport == "BOINC"):
            msg = "writeCntlFile: invalid JOB_TRANSPORT: %s Choose one of 'SSH', 'BOINC'." % job_transport
            self.exit(msg)
        input += "JOB_TRANSPORT = '%s'\n" % job_transport
            
        re_type = self.keywords.get('RE_TYPE')
        if re_type is None:
            msg = "writeCntlFile: RE_TYPE is not specified"
            self.exit(msg)
        if not (re_type == 'TEMPT' or re_type == 'BEDAMTEMPT'):
            msg = "writeCntlFile: invalid RE_TYPE. Choose one of 'TEMPT', 'BEDAMTEMPT'."
            self.exit(msg)
        input += "RE_TYPE = '%s'\n" % re_type

        engine = self.keywords.get('ENGINE')
        if engine is None:
            engine = "IMPACT"
        input += "ENGINE = '%s'\n" % engine

        input += "ENGINE_INPUT_BASENAME = '%s'\n" % self.jobname

        if job_transport == 'SSH':
            multiarch = self.keywords.get('MULTIARCH')
            if multiarch == 'YES':
                exec_directory =  self.keywords.get('EXEC_DIRECTORY')
                if exec_directory is None:
                    msg = "writeCntlFile: multiarch ASyncRE requires EXEC_DIRECTORY"
                    self.exit(msg)
                input += "EXEC_DIRECTORY = '%s'\n" % exec_directory

        input += "RE_SETUP = 'YES'\n"

        extfiles = self.keywords.get('ENGINE_INPUT_EXTFILES')
        required_files = "runimpact"
        if extfiles is None:
            extfiles = required_files
        else:
            extfiles += ",%s" % required_files
        restart_file = "%s_0.rst" % self.jobname
        input_file = "%s.inp" % self.jobname
        extfiles = extfiles + ",%s,%s" % (restart_file,input_file)
        if re_type == 'BEDAMTEMPT':
            extfiles += ",%s,%s,%s" % (self.recidxfile,self.ligidxfile,self.restraint_file)
        input += "ENGINE_INPUT_EXTFILES = '%s'\n" % extfiles
        
        temperatures = self.keywords.get('TEMPERATURES')
        if temperatures is not None:
            input += "TEMPERATURES = '%s'\n" % temperatures
        
        lambdas = self.keywords.get('LAMBDAS')
        if lambdas is not None:
            input += "LAMBDAS = '%s'\n" % lambdas

        wall_time = self.keywords.get('WALL_TIME')
        if wall_time is not None:
            input += "WALL_TIME = %d\n" % int(wall_time)

        replica_run_time = self.keywords.get('REPLICA_RUN_TIME')
        if replica_run_time is not None:
            input += "REPLICA_RUN_TIME = %d\n" % int(replica_run_time)

        cycle_time = self.keywords.get('CYCLE_TIME')
        if cycle_time is not None:
            input += "CYCLE_TIME = %d\n" % int(cycle_time)
            
        if job_transport == 'SSH':
            input += "NODEFILE = 'nodefile'\n"
            total_cores = self.keywords.get('TOTAL_CORES')
            if total_cores is None:
                msg = "writeCntlFile: TOTAL_CORES is required"
                self.exit(msg)
            input += "TOTAL_CORES = %d\n" % int(total_cores) 
            subjob_cores = self.keywords.get('SUBJOB_CORES')
            if subjob_cores is not None:
                input += "SUBJOB_CORES = %d\n" % int(subjob_cores)

        if job_transport == 'BOINC':
            total_cores = self.keywords.get('TOTAL_CORES')
            if total_cores is None:
                msg = "writeCntlFile: TOTAL_CORES is required"
                self.exit(msg)
            input += "TOTAL_CORES = %d\n" % int(total_cores) 
            subjob_cores = self.keywords.get('SUBJOB_CORES')
            if subjob_cores is not None:
                input += "SUBJOB_CORES = %d\n" % int(subjob_cores)

        if job_transport == 'BOINC':
            boinc_projectdir = self.keywords.get('BOINC_PROJECTDIR')
            if boinc_projectdir is None:
                msg = "writeCntlFile: BOINC_PROJECTDIR is required"
                self.exit(msg)
            boinc_database = self.keywords.get('BOINC_DATABASE')
            if boinc_database is None:
                msg = "writeCntlFile: BOINC_DATABASE is required"
                self.exit(msg)
            boinc_database_user = self.keywords.get('BOINC_DATABASE_USER')
            if boinc_database_user is None:
                msg = "writeCntlFile: BOINC_DATABASE_USER is required"
                self.exit(msg)
            boinc_database_password = self.keywords.get('BOINC_DATABASE_PASSWORD')
            if boinc_database_password is None:
                msg = "writeCntlFile: BOINC_DATABASE_PASSWORD is required"
                self.exit(msg)
            input += "BOINC_PROJECTDIR = '%s'\n" % boinc_projectdir
            input += "BOINC_DATABASE = '%s'\n" % boinc_database
            input += "BOINC_DATABASE_USER = '%s'\n" % boinc_database_user
            input += "BOINC_DATABASE_PASSWORD = '%s'\n" % boinc_database_password

        subjobs_buffer_size = self.keywords.get('SUBJOBS_BUFFER_SIZE')
        if subjobs_buffer_size is not None:
            input += "SUBJOBS_BUFFER_SIZE = '%f'\n" % float(subjobs_buffer_size)

        verbose = self.keywords.get('VERBOSE')
        if verbose is not None:
            input += "VERBOSE = '%s'\n" % verbose
        
        cntlfile = "%s_asyncre.cntl" % self.jobname
        f = open(cntlfile, "w")
        f.write(input)
        f.close

#
# writes the Impact input file for minimization/thermalization
#
    def  writeThermInputFile(self):
        if self.recidxfile is None:
            msg = "writeThermInputFile: Internal error: receptor structure file not found"
            self.exit(msg)
        if not os.path.exists(self.recidxfile):
            msg = 'File does not exist: %s' % self.recidxfile
            self.exit(msg)
        if self.ligidxfile is None:
            msg = "writeThermInputFile: Internal error: receptor structure file not found"
            self.exit(msg)
        if not os.path.exists(self.ligidxfile):
            msg = 'File does not exist: %s' % self.ligidxfile
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
        out_restart_file =    self.jobname + '_0' + '.rst'
        self.mintherm_out_restart_file = out_restart_file 

        out_rcpt_structure_file =  self.jobname + '_rcpt_mintherm' + '.dms'
        out_lig_structure_file =  self.jobname + '_lig_mintherm' + '.dms'

        input = self.input_mintherm.format(
            out_file = impact_output_file, title = impact_jobtitle,
            dms_rcpt_in = self.recidxfile, dms_lig_in = self.ligidxfile,
            cmkf = kfcm, cmdist0 = d0cm, cmtol = tolcm, cmrestraints_file = self.restraint_file,
            temperature = temperature,
            rst_file_out = out_restart_file,
            dms_rcpt_out =  out_rcpt_structure_file,
            dms_lig_out = out_lig_structure_file)

        f = open(impact_input_file, "w")
        f.write(input)
        f.close


#
# writes the Impact input file for AsyncRE production
#
    def  writeRemdInputFile(self):
        if self.recidxfile is None:
            msg = "writeThermInputFile: Internal error: receptor structure file not found"
            self.exit(msg)
        if not os.path.exists(self.recidxfile):
            msg = 'File does not exist: %s' % self.recidxfile
            self.exit(msg)
        if self.ligidxfile is None:
            msg = "writeThermInputFile: Internal error: receptor structure file not found"
            self.exit(msg)
        if not os.path.exists(self.ligidxfile):
            msg = 'File does not exist: %s' % self.ligidxfile
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

        rest_kf = self.keywords.get('REST_RECEPTOR_KF')
        if not rest_kf:
            rest_kf = '0.6'

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

        input = self.input_remd.format(
            job_name = self.jobname,
            dms_rcpt_in = self.recidxfile, dms_lig_in = self.ligidxfile,
            umax0 = umax, rest_kf = rest_kf,
            cmkf = kfcm, cmdist0 = d0cm, cmtol = tolcm, cmrestraints_file = self.restraint_file,
            nmd_eq = nmd_eq, 
            nmd = nmd_prod, nprnt = nprnt)

        impact_input_file = self.jobname + ".inp"

        f = open(impact_input_file, "w")
        f.write(input)
        f.close

#
# writes the 'runimpact' script for interactive and asyncre use
#
    def  writeRunimpactFile(self):
        
        job_transport = self.keywords.get('JOB_TRANSPORT')
        if job_transport is None:
            msg = "writeRunimpactFile: JOB_TRANSPORT is not specified"
            self.exit(msg)

        if job_transport == 'SSH':
            multiarch = self.keywords.get('MULTIARCH')
            if multiarch is not None and multiarch == 'YES':
                input = self.input_runimpact_multiarch
            else:
                subjob_cores = self.keywords.get('SUBJOB_CORES')
                if subjob_cores is None:
                    subjob_cores = 1
                exec_directory =  self.keywords.get('EXEC_DIRECTORY')
                if exec_directory is None:
                    msg = "writeRunimpactFile: EXEC_DIRECTORY is required"
                    self.exit(msg)
                input = self.input_runimpact_standard.format(
                    acd_impact_home = exec_directory,
                    subjob_cores = subjob_cores)
            f = open('runimpact', "w")
            f.write(input)
            f.close
            os.chmod('runimpact',0744)

        #this is the interactive version
        input = self.input_runimpact_interactive
        f = open('runimpact_i', "w")
        f.write(input)
        f.close
        os.chmod('runimpact_i',0744)

#
# writes sample submission scripts for slurm (stampede) and PBS
#
    def writeQueueFiles(self):
        input = self.input_slurm.format(job_name = self.jobname)
        slurm_file = self.jobname + '.slurm'
        f = open(slurm_file, "w")
        f.write(input)
        f.close

        input = self.input_qsub.format(job_name = self.jobname)
        qsub_file = self.jobname + '.qsub'
        f = open(qsub_file, "w")
        f.write(input)
        f.close

##################### MAIN CODE ##########################
if __name__ == '__main__':

    # Setup the logger
    logger = schrodinger.utils.log.get_output_logger("bedam_prep")

    # Parse arguments:
    usage = "%prog [options] <inputfile>"
    parser = cmdline.SingleDashOptionParser(usage)
    (options, args) = parser.parse_args(sys.argv[1:])
    
    if len(args) != 1:
        parser.error("Please specify ONE input file")
    
    commandFile = args[0]

    print ""
    print "===================================="
    print "       BEDAM Job Preparation        "
    print "===================================="
    print ""
    print "SCHRODINGER: " + os.environ['SCHRODINGER']
    print "Started at: " + str(time.asctime())
    print "Input file:", commandFile
    print ""
    sys.stdout.flush()
    
    print "Reading options"
    bedam = bedam_job_asyncre(commandFile, options)

    print "Set put templates for input files ..."
    bedam.setupTemplatesASyncRE()
    print "Analyzing structure files ..."
    bedam.getDesmondDMSFiles()
    print "Writing BEDAM restraint file ..."
    bedam.writeRestraintFile()
    print "Adding atomic restraints to receptor file ..."
    bedam.writeRecStructureFile()
    print "Writing job input files ..."
    bedam.writeCntlFile()
    bedam.writeThermInputFile()    
    bedam.writeRemdInputFile()
    bedam.writeRunimpactFile()
    bedam.writeQueueFiles()

    print
    print "Job preparation complete"
    print ""
    print "To run the minimization/thermalization calculation do:"
    exec_directory =  bedam.keywords.get('EXEC_DIRECTORY')
    if exec_directory is not None:
        print "export IMPACTHOME=" + exec_directory
    else:
        print "export IMPACTHOME=" + "<path_to_academic_impact_directory>"
    print "export OMP_NUM_THREADS=" + "<number_of_CPU_cores>"
    print "./runimpact_i %s_mintherm.inp" % bedam.jobname
    print ""
    print "When completed run the production calculation with:"
    print "<path_to_asyncre>/bedamtempt_async_re.py %s_asyncre.cntl" % bedam.jobname
    print ""
