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

#
#  write the dihedral list file for flattening simulation 
#
    def writeDihedralFlatFile(self):
        # check that structure files with internal indexes have been generated
        if self.recidxfile is None or self.ligidxfile is None:
            msg = "writeRestraintFile: Internal error: structure files not found"
            self.exit(msg)

        rcp_dihe_str =  self.keywords.get('FLAT_DIHE_RECSQL')
        lig_dihe_str =  self.keywords.get('FLAT_DIHE_LIGSQL')

        if rcp_dihe_str or lig_dihe_str :
           dihe_flatfile = self.keywords.get('FLAT_DIHE_FILE')
           f = open(dihe_flatfile,"w")

        if not rcp_dihe_str:
            msg = "bedam_prep: No atom selection specified for dihedral flattening of receptor"
            print msg
        else : 
            atompairs = rcp_dihe_str.split(',')
            npairs = len(atompairs)
            sec_atoms=[]
            thd_atoms=[]
            for atompair in atompairs :
                tempatoms = atompair.split(' ')
                sec_atoms.append(int(tempatoms[0])-1)
                thd_atoms.append(int(tempatoms[1])-1)
            con = lite.connect(self.recidxfile)
            ii_atoms_idx=[] 
            with con: 
                cur = con.cursor()
                cur.execute("SELECT i_i_internal_atom_index FROM particle;")
                rows=cur.fetchall()
                for row in rows:
                    ii_atoms_idx.append(row[0])

                flat_dihlist=[]
                for i in range(0,len(sec_atoms)) :
                    cur.execute("SELECT dihedral_trig_term.p0,dihedral_trig_term.p1,dihedral_trig_term.p2,dihedral_trig_term.p3,dihedral_trig_term.param,pair_12_6_es_term.param FROM dihedral_trig_term LEFT OUTER JOIN pair_12_6_es_term ON (dihedral_trig_term.p0=pair_12_6_es_term.p0 AND dihedral_trig_term.p3=pair_12_6_es_term.p1) OR (dihedral_trig_term.p3=pair_12_6_es_term.p0 AND dihedral_trig_term.p0=pair_12_6_es_term.p1) where (dihedral_trig_term.p1=? AND dihedral_trig_term.p2=?);", (sec_atoms[i],thd_atoms[i]))
                    rows=cur.fetchall()
                    for row in rows: 
                       if row[5] is not None :
                          p0t3=[]
                          p0t3.append(row[0])
                          p0t3.append(row[1])
                          p0t3.append(row[2])
                          p0t3.append(row[3])
                          flat_dihlist.append(p0t3)
                    vn=[]
                    phase=[1,1,-1,1,-1,1,-1]
                    for i in range(0,7) :
                        vn.append(0.0)

                for j in range(0, len(flat_dihlist)) :
                    cur.execute("SELECT p0,p1,p2,p3,phi0,fc0,fc1,fc2,fc3,fc4,fc5,fc6 FROM dihedral_trig where (p0=? AND p1=? AND p2=? AND p3=?)", (flat_dihlist[j][0],flat_dihlist[j][1],flat_dihlist[j][2],flat_dihlist[j][3]))
	            rows=cur.fetchall()
                    for row in rows:
                        sum = 0.0
                        for i in range(1,7) :
                            vn[i]=row[i+5]*phase[i]
                            sum += row[i+5]*phase[i]
                        vn[0]=row[5]-sum
                        i_int=int(ii_atoms_idx[row[0]])
                        j_int=int(ii_atoms_idx[row[1]])
                        k_int=int(ii_atoms_idx[row[2]])
                        l_int=int(ii_atoms_idx[row[3]])
                        if (abs(vn[0])>0.0) :
                           f.write("%d %d %d %d %d %f %d %d\n" %(1,i_int,j_int,k_int,l_int,vn[0],phase[0],0))
                        for i in range(1,7) :
                            if abs(vn[i]) > 0.0 : 
                               f.write("%d %d %d %d %d %f %d %d\n" %(1,i_int,j_int,k_int,l_int,vn[i],phase[i],i))
            con.close()

        if not lig_dihe_str:
            msg = "bedam_prep: No atom selection specified for dihedral flattening of ligand"
            print msg
        else :
            atompairs = lig_dihe_str.split(',')
            npairs = len(atompairs)
            sec_atoms=[]
            thd_atoms=[]
            for atompair in atompairs :
                tempatoms = atompair.split(' ')
                sec_atoms.append(int(tempatoms[0])-1)
                thd_atoms.append(int(tempatoms[1])-1)
            con = lite.connect(self.ligidxfile)
            ii_atoms_idx=[]
            with con:
                cur = con.cursor()
                cur.execute("SELECT i_i_internal_atom_index FROM particle;")
                rows=cur.fetchall()
                for row in rows:
                    ii_atoms_idx.append(row[0])

                flat_dihlist=[]
                for i in range(0,len(sec_atoms)) :
                    cur.execute("SELECT dihedral_trig_term.p0,dihedral_trig_term.p1,dihedral_trig_term.p2,dihedral_trig_term.p3,dihedral_trig_term.param,pair_12_6_es_term.param FROM dihedral_trig_term LEFT OUTER JOIN pair_12_6_es_term ON (dihedral_trig_term.p0=pair_12_6_es_term.p0 AND dihedral_trig_term.p3=pair_12_6_es_term.p1) OR (dihedral_trig_term.p3=pair_12_6_es_term.p0 AND dihedral_trig_term.p0=pair_12_6_es_term.p1) where (dihedral_trig_term.p1=? AND dihedral_trig_term.p2=?);", (sec_atoms[i],thd_atoms[i]))
                    rows=cur.fetchall()
                    for row in rows:
                       if row[5] is not None :
                          p0t3=[]
                          p0t3.append(row[0])
                          p0t3.append(row[1])
                          p0t3.append(row[2])
                          p0t3.append(row[3])
                          flat_dihlist.append(p0t3)
                    vn=[]
                    phase=[1,1,-1,1,-1,1,-1]
                    for i in range(0,7) :
                        vn.append(0.0)

                for j in range(0, len(flat_dihlist)) :
                    cur.execute("SELECT p0,p1,p2,p3,phi0,fc0,fc1,fc2,fc3,fc4,fc5,fc6 FROM dihedral_trig where (p0=? AND p1=? AND p2=? AND p3=?)", (flat_dihlist[j][0],flat_dihlist[j][1],flat_dihlist[j][2],flat_dihlist[j][3]))
                    rows=cur.fetchall()
                    for row in rows:
                        sum = 0.0
                        for i in range(1,7) :
                            vn[i]=row[i+5]*phase[i]
                            sum += row[i+5]*phase[i]
                        vn[0]=row[5]-sum
                        i_int=int(ii_atoms_idx[row[0]])
                        j_int=int(ii_atoms_idx[row[1]])
                        k_int=int(ii_atoms_idx[row[2]])
                        l_int=int(ii_atoms_idx[row[3]])
                        if (abs(vn[0])>0.0) :
                           f.write("%d %d %d %d %d %f %d %d\n" %(2,i_int,j_int,k_int,l_int,vn[0],phase[0],0))
                        for i in range(1,7) :
                            if abs(vn[i]) > 0.0 :
                               f.write("%d %d %d %d %d %f %d %d\n" %(2,i_int,j_int,k_int,l_int,vn[i],phase[i],i))
            con.close()
        f.close()
                
#
#  write the atom list file for nonbonded flattening simulation 
#
    def writeNonbFlatFile(self):
        # check that structure files with internal indexes have been generated
        if self.recidxfile is None or self.ligidxfile is None:
            msg = "writeRestraintFile: Internal error: structure files not found"
            self.exit(msg)

        rcp_nonb_str =  self.keywords.get('FLAT_NONB_RECSQL')
        lig_nonb_str =  self.keywords.get('FLAT_NONB_LIGSQL')

        if rcp_nonb_str or lig_nonb_str :
           nonb_flatfile = self.keywords.get('FLAT_NONB_FILE')
           f = open(nonb_flatfile,"w")

        if not rcp_nonb_str:
            msg = "bedam_prep: No atom selection specified for nonbonded flattening of receptor"
            print msg
        else :
            atoms = rcp_nonb_str.split(',')
            natoms = len(atoms)
            nonb_atoms=[]
            for atom in atoms :
                nonb_atoms.append(int(atom)-1)
            con = lite.connect(self.recidxfile)
            ii_atoms_idx=[]
            with con:
                cur = con.cursor()
                cur.execute("SELECT i_i_internal_atom_index FROM particle;")
                rows=cur.fetchall()
                for row in rows:
                    ii_atoms_idx.append(row[0])
                for i in range(0,natoms) :
                    i_int=int(ii_atoms_idx[nonb_atoms[i]])
                    f.write("%d %d\n" %(1,i_int))
            con.close()

        if not lig_nonb_str:
            msg = "bedam_prep: No atom selection specified for nonbonded flattening of ligand"
            print msg
        else :
            atoms = lig_nonb_str.split(',')
            natoms = len(atoms)
            nonb_atoms=[]
            for atom in atoms :
                nonb_atoms.append(int(atom)-1)
            con = lite.connect(self.ligidxfile)
            ii_atoms_idx=[]
            with con:
                cur = con.cursor()
                cur.execute("SELECT i_i_internal_atom_index FROM particle;")
                rows=cur.fetchall()
                for row in rows:
                    ii_atoms_idx.append(row[0])
                for i in range(0,natoms) :
                    i_int=int(ii_atoms_idx[nonb_atoms[i]])
                    f.write("%d %d\n" %(2,i_int))
            con.close()
        f.close() 
    
##################### MAIN CODE ##########################
if __name__ == '__main__':

    # Setup the logger
    logger = schrodinger.utils.log.get_output_logger("bedam_prep_acflat")

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
    print "Create dihedral list for flattening  ..."
    bedam.writeDihedralFlatFile()
    print "Create nonbonded list for flattening  ..."
    bedam.writeNonbFlatFile()
    print "Writing Minimization/Thermalization input file ..."
    bedam.writeThermInputFile()
    print "Writing Remd Production input file ..."
    bedam.writeRemdInputFile()
    print "Writing Trajectory conversion input file ..."
    bedam.writeReadTrajInputFile()

    print
    print "Job preparation complete"

#EOF
