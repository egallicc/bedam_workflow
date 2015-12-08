# bedam python class
__doc__="""
$Revision: 0.1 $

A class to prepare BEDAM jobs using academic IMPACT

"""
# Contributors: Junchao Xia

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
from  bedam import bedam_job  


class bedam_prep_ac (bedam_job):
    """
    Class to set up remd calculations
    """
    def __init__(self, command_file, options):
        bedam_job.__init__(self, command_file, options)
        self.setupTemplates()	

#
# Impact input file templates for academic
#
    def setupTemplates(self):
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
        self.input_agbnp2 =  """
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
    mmechanics consolv agbnp2
    weight constraints buffer 25.000000
  quit
  read parm file -
"paramstd.dat" -
  noprint
  energy parm dielectric 1 nodist -
   listupdate 10 -
    cutoff 12
  energy rescutoff byatom all
  zonecons auto
  energy constraints bonds hydrogens
QUIT

MINIMIZE
  input cntl mxcyc 0  rmscut 0.05 deltae 1.0e-05
  conjugate dx0 0.05 dxm 1.0
  run
  write sql name species1 file "%s"
  write sql name species2 file "%s"
QUIT

END
"""
        self.input_agbnp2_sep =  """
write file -
"%s" -
      title -
"%s" *

CREATE
  build primary name species1 type auto read sqldb file -
"%s"
QUIT

SETMODEL
  setpotential
    mmechanics consolv agbnp2
    weight constraints buffer 25.000000
  quit
  read parm file -
"paramstd.dat" -
  noprint
  energy parm dielectric 1 nodist -
   listupdate 10 -
    cutoff 12
  energy rescutoff byatom all
  zonecons auto
  energy constraints bonds hydrogens
QUIT

MINIMIZE
  input cntl mxcyc 0  rmscut 0.05 deltae 1.0e-05
  conjugate dx0 0.05 dxm 1.0
  run
  write sql name species1 file "%s"
QUIT

END
"""
        self.input_idx =  """
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
  write sql name species1 file "%s"
  write sql name species2 file "%s"
QUIT

END
"""
        self.input_mintherm = """
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
  write sql name species1 file "%s"
  write sql name species2 file "%s"
QUIT

END
"""
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
  build primary name species1 type auto read sqldb file -
"%s"
  build primary name species2 type auto read sqldb file -
"%s"
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
# runs desmond utils from commercial schrodinger package on the input mae files to get versions with internal 
# atom indexes
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
        rcpt_cms_file = self.jobname + '_rcpt-out.cms'
        lig_cms_file = self.jobname + '_lig-out.cms'
        f = open(desmond_builder_file, 'w')
        input =  self.input_cms
        f.write(input)
        f.close()
        source_cmd = '. ' + com_shrod_source
        rcpt_cmd = '$SCHRODINGER/utilities/multisim' + ' -JOBNAME ' + self.jobname + ' -m ' + desmond_builder_file + ' ' + receptor_file + ' -o ' + rcpt_cms_file + ' -HOST localhost -maxjob 1 -WAIT'
        lig_cmd =  '$SCHRODINGER/utilities/multisim' + ' -JOBNAME ' + self.jobname + ' -m ' + desmond_builder_file + ' ' + ligand_file + ' -o ' + lig_cms_file + ' -HOST localhost -maxjob 1 -WAIT' 
        # os.system(source_cmd)
        # os.system(rcpt_cmd)
        # os.system(lig_cmd)
        cms_cmd = source_cmd + ";" + rcpt_cmd + ";" + lig_cmd 
        os.system(cms_cmd)

        print "Convert cms files to dms files"
        rcpt_dms_file = self.jobname + '_rcpt-out.dms'
        lig_dms_file = self.jobname + '_lig-out.dms'
	rcpt_cmd = '$SCHRODINGER/run -FROM desmond mae2dms ' + rcpt_cms_file + ' ' + rcpt_dms_file
        lig_cmd = '$SCHRODINGER/run -FROM desmond mae2dms ' + lig_cms_file + ' ' + lig_dms_file
        # os.system(source_cmd)
        # os.system(rcpt_cmd)
        # os.system(lig_cmd)
        dms_cmd = source_cmd + ";" + rcpt_cmd + ";" + lig_cmd
        os.system(dms_cmd)


        print "Add agbnp parameters into dms files"

        self.agbnp_script = self.keywords.get('AGBNP_SCRIPT')
        if self.agbnp_script is None:
            self.exit('AGBNP_SCRIPT needs to be specified to add agbnp parameter.')

        agbnp_receptor_file =   self.jobname + '_rcpt_agbnp' + '.dms'
        copy_cmd = "cp " + rcpt_dms_file + " " + agbnp_receptor_file
        agbnp_cmd =  "$SCHRODINGER/run " + self.agbnp_script + "/add_agbnp2.py " + agbnp_receptor_file
        agbnp_cmd = source_cmd + ";" + copy_cmd + ";" + agbnp_cmd
        os.system(agbnp_cmd)

        agbnp_ligand_file =     self.jobname + '_lig_agbnp' + '.dms'
        copy_cmd = "cp " + lig_dms_file + " " + agbnp_ligand_file
        agbnp_cmd =  "$SCHRODINGER/run " + self.agbnp_script + "/add_agbnp2.py " + agbnp_ligand_file
        agbnp_cmd = source_cmd + ";" + copy_cmd + ";" + agbnp_cmd
        os.system(agbnp_cmd)

        # Add agbnp parameters using commerical impact 
        #print "add agbnp parameters into dms files"
        #agbnp_shrod_source =  self.keywords.get('AGBNP_IMPACT_EVN')
        #if not agbnp_shrod_source:
        #    msg = "bedam_prep: No commercial IMPACT source file for agbnp specified in the input file"
        #    self.exit(msg) 
        # 
        #agbnp_input_file =   self.jobname + '_rcpt_agbnp' + '.inp'
        #agbnp_output_file =  self.jobname + '_rcpt_agbnp' + '.out'
        #agbnp_jobtitle =     self.jobname + '_rcpt_agbnp'
        #agbnp_receptor_file =   self.jobname + '_rcpt_agbnp' + '.dms'

        #f = open(agbnp_input_file, 'w')
        #input =  self.input_agbnp2_sep % (agbnp_output_file, agbnp_jobtitle, rcpt_dms_file, agbnp_receptor_file)
        #f.write(input)
        #f.close()
        #agbnp_log_file =  self.jobname + '_rcpt_agbnp' + '.log'
        #source_cmd = '. ' + agbnp_shrod_source 
        # agbnp_cmd =  "$SCHRODINGER/impact -i " + agbnp_input_file  + " -WAIT -LOCAL "
        #agbnp_cmd =  "$IMPACT_EXEC/main1m " + agbnp_input_file  + " > " + agbnp_log_file + " 2>&1 "
        #agbnp_cmd = source_cmd + ";" + agbnp_cmd
        #os.system(agbnp_cmd)

        #agbnp_input_file =   self.jobname + '_lig_agbnp' + '.inp'
        #agbnp_output_file =  self.jobname + '_lig_agbnp' + '.out'
        #agbnp_jobtitle =     self.jobname + '_lig_agbnp'
        #agbnp_ligand_file =     self.jobname + '_lig_agbnp' + '.dms'
        #f = open(agbnp_input_file, 'w')
        #input =  self.input_agbnp2_sep % (agbnp_output_file, agbnp_jobtitle,lig_dms_file, agbnp_ligand_file )
        #f.write(input)
        #f.close()
        #agbnp_log_file =  self.jobname + '_lig_agbnp' + '.log'
        #source_cmd = '. ' + agbnp_shrod_source 
        ## agbnp_cmd =  "$SCHRODINGER/impact -i " + agbnp_input_file  + " -WAIT -LOCAL "
        #agbnp_cmd =  "$IMPACT_EXEC/main1m " + agbnp_input_file  + " > " + agbnp_log_file + " 2>&1 "
        #agbnp_cmd = source_cmd + ";" + agbnp_cmd
        #os.system(agbnp_cmd)

        #if not os.path.exists(agbnp_receptor_file) or not os.path.exists(agbnp_ligand_file):
        #    print " failed"
        #    msg = "Impact job to generate agbnp dms files failed"
        #    self.exit(msg)

        print "add internal atom indices into dms files"
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
        input =  self.input_idx % (impact_output_file, impact_jobtitle, agbnp_receptor_file, agbnp_ligand_file, out_receptor_file, out_ligand_file )

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
        receptor_sql =  self.keywords.get('REST_LIGAND_CMRECSQL')
        if not receptor_sql:
            msg = "bedam_prep: No sql selection specified for receptor site center"
            self.exit(msg)

        ligand_sql =  self.keywords.get('REST_LIGAND_CMLIGSQL')
        if not ligand_sql:
            ligand_sql = '( all)'
        recpt_cmd = "SELECT id,x,y,z,i_i_internal_atom_index FROM particle WHERE " + receptor_sql
 
        con = lite.connect(self.recidxfile)
        with con:
            cur = con.cursor()  
            cur.execute(recpt_cmd)
            rows = cur.fetchall()
            if not rows:
                msg = "bedam_prep: No atoms return from sql selection for receptor site center"
                self.exit(msg) 
            rec_atoms = []
            cmrx = cmry = cmrz = 0. 
            for row in rows:
                cmrx += row[1]
                cmry += row[2]
                cmrz += row[3]
                rec_atoms.append(row[4])
            n=len(rows)
            cmrx /= float(n)
            cmry /= float(n)
            cmrz /= float(n)

        lig_cmd = "SELECT id,x,y,z,i_i_internal_atom_index FROM particle WHERE " + ligand_sql 
        con = lite.connect(self.ligidxfile)
        with con:
            cur = con.cursor()  
            cur.execute(lig_cmd)
            rows = cur.fetchall()
            if not rows:
                msg = "bedam_prep: No atoms return from sql selection for ligand site center"
                self.exit(msg) 
            lig_atoms = []
            cmlx = cmly = cmlz = 0. 
            for row in rows:
                cmlx += row[1]
                cmly += row[2]
                cmlz += row[3]
                lig_atoms.append(row[4])
            n=len(rows)
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
# writes receptor dms file with restraints
#
    def writeRecStructureFile(self):
        if self.recidxfile is None :
            msg = "writeRecStructureFile: Internal error: Structure file not found"
            self.exit(msg)            
        if not os.path.exists(self.recidxfile):
            msg = 'File does not exist: %s' % self.recidxfile
            self.exit(msg)

        rest_sql =  self.keywords.get('REST_RECEPTOR_SQL')
        if rest_sql is not None: 
            con = lite.connect(self.recidxfile)
            recpt_cmd = "PRAGMA table_info(particle)";
            columnExists = False; 
            with con:
                cur = con.cursor()  
                cur.execute(recpt_cmd)
                rows = cur.fetchall()
                for row in rows:
                    if row[1] == "grp_buffer" : 
                        columnExists = True
                if not columnExists : 
                    cur.execute("ALTER TABLE particle ADD COLUMN grp_buffer int DEFAULT(0);")
                recpt_cmd = "SELECT id FROM particle WHERE " + rest_sql
                cur.execute(recpt_cmd)
                atoms = cur.fetchall()      
                for iat in atoms:
                    cur.execute("UPDATE particle SET grp_buffer=? WHERE Id=?", (2, iat[0]))
 
        froz_sql =  self.keywords.get('FROZ_RECEPTOR_SQL')
        if froz_sql is not None: 
            con = lite.connect(self.recidxfile)
            recpt_cmd = "PRAGMA table_info(particle)";
            columnExists = False; 
            with con:
                cur = con.cursor()  
                cur.execute(recpt_cmd)
                rows = cur.fetchall()
                for row in rows:
                    if row[1] == "grp_frozen" : 
                        columnExists = True
                if not columnExists : 
                    cur.execute("ALTER TABLE particle ADD COLUMN grp_frozen int DEFAULT(0);")
                recpt_cmd = "SELECT id FROM particle WHERE " + froz_sql
                cur.execute(recpt_cmd)
                atoms = cur.fetchall()      
                for iat in atoms:
                    cur.execute("UPDATE particle SET grp_frozen=? WHERE Id=?", (1, iat[0]))

        receptor_file_restr =   self.jobname + '_rcpt_restr' + '.dms'
        cpcmd = "cp " +  self.recidxfile + "  " + receptor_file_restr 
        os.system(cpcmd)
        self.receptor_file_restr = receptor_file_restr
        self.ligand_file_restr = self.ligidxfile

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
        if self.ligand_file_restr is None:
            msg = "writeThermInputFile: No ligand file specified in the input file"
            self.exit(msg)
        if not os.path.exists(self.ligand_file_restr):
            msg = 'File does not exist: %s' % self.ligand_file_restr
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
        self.mintherm_out_restart_file = out_restart_file 
        out_rcpt_structure_file =  self.jobname + '_rcpt_mintherm' + '.dms'
        out_lig_structure_file =  self.jobname + '_lig_mintherm' + '.dms'
        input = self.input_mintherm % (impact_output_file, impact_jobtitle, self.receptor_file_restr,self.ligand_file_restr, kfcm, d0cm, tolcm, self.restraint_file, temperature, out_restart_file, out_rcpt_structure_file,out_lig_structure_file)
        f = open(impact_input_file, "w")
        f.write(input)
        f.close()

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
        if self.ligand_file_restr is None:
            msg = "writeRemdInputFile: No ligand file specified in the input file"
            self.exit(msg)
        if not os.path.exists(self.ligand_file_restr):
            msg = 'File does not exist: %s' % self.ligand_file_restr
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
        out_rcpt_structure_file =  self.jobname + '_rcpt_remd' + '.dms'
        out_lig_structure_file =  self.jobname + '_lig_remd' + '.dms'
        if not restart_file:
           input = self.input_remd % (impact_output_file, impact_jobtitle, 
                                      self.receptor_file_restr , self.ligand_file_restr,
                                      umax, bias_flag,
                                      rest_kf, nlambdas, lambda_list,
                                      bias_parameters,
                                      kfcm, d0cm, tolcm, self.restraint_file, 
                                      self.mintherm_out_restart_file,
                                      nmd_eq, nlambdas, nmd_eq, temperature, temperature, nprnt,
                                      nmd_prod, nlambdas, rsvr_command,
                                      temperature, temperature, nprnt,
                                      ntrj, out_trajectory_file, out_restart_file, 
                                      out_rcpt_structure_file, out_lig_structure_file)
        else:
           input = self.input_remd_restart % (impact_output_file, impact_jobtitle, 
                                              self.receptor_file_restr , self.ligand_file_restr,
                                              umax, bias_flag,
                                              rest_kf, nlambdas, lambda_list,
                                              bias_parameters,
                                              kfcm, d0cm, tolcm, self.restraint_file, 
                                              nmd_prod, nlambdas, rsvr_command,
                                              temperature, temperature, nprnt,
                                              restart_file,
                                              ntrj, out_trajectory_file, 
                                              out_restart_file, out_rcp_structure_file, out_lig_structure_file)
                
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
        if self.ligand_file_restr is None:
            msg = "writeReadTrajInputFile: No ligand file specified in the input file"
            self.exit(msg)
        if not os.path.exists(self.ligand_file_restr):
            msg = 'File does not exist: %s' % self.ligand_file_restr
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
                                   self.receptor_file, self.ligand_file, 
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
                                   self.receptor_file_restr , self.ligand_file_restr, 
                                  umax,
                                   rest_kf, nlambdas, lambda_list, 
                                   nlambdas,nlambdas-1)
        f = open(impact_input_file, "w")
        f.write(input)
        f.close()

    
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
    bedam = bedam_prep_ac(commandFile, options)
    print "Reading the templates for input files ..."
    bedam.setupTemplates()
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
