# python script
__doc__="""
$Revision: 0.1 $

BEDAM analysis script

"""
# Contributors:  Ahmet Mentes and  Junchao Xia

import os, sys, time
from operator import itemgetter, attrgetter
from schrodinger import structure, structureutil
from schrodinger.application import inputconfig
from schrodinger.utils import cmdline
import schrodinger.utils.log

from math import *
from numpy import * # numerical array library

from bedam import bedam_job


class bedam_analyze_acflat (bedam_job):
    """
    Class to set up remd calculations
    """
    def __init__(self, command_file, options):
        bedam_job.__init__(self, command_file, options)


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
                lmb = datai[i][nf-4]
                u = datai[i][nf-3]
                be.append(u)
                lv.append(lmb)
                if not lmb in self.binding_energies_bylambda.keys():
                    self.binding_energies_bylambda[lmb] = []
                self.binding_energies_bylambda[lmb].append(u)
            self.binding_energies_bysim.append(be)
            self.lambda_bysim.append(lv)

    def getFlatteningEnergies(self):
        """
        extract flattening energies from .out files
        """
        lambda_string = self.keywords.get('LAMBDAS')
        if not lambda_string:
            msg = "getFlatteningEnergies: No lambda schedule specified"
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

        self.flattening_energies_bysim = []
        self.lambda_bysim = []
        self.flattening_energies_bylambda = {}
        for i in range(0,nlambdas):
            impact_output_file =  self.jobname + '_remd' '_%d' % i + '.out'
            datai = self.getImpactData(impact_output_file)
            fe = []
            lv = []
            ndata = len(datai)
            nf = len(datai[0])
            for i in range(nreject,ndata):
                lmb = datai[i][nf-4]
                eflat = datai[i][nf-2]
                fe.append(eflat)
                lv.append(lmb)
                if not lmb in self.flattening_energies_bylambda.keys():
                    self.flattening_energies_bylambda[lmb] = []
                self.flattening_energies_bylambda[lmb].append(eflat)
            self.flattening_energies_bysim.append(fe)
            self.lambda_bysim.append(lv)

    def getFlatteningEnergies2(self):
        """
        extract flattening energies from .out files
        """
        lambda_string = self.keywords.get('LAMBDAS')
        if not lambda_string:
            msg = "getFlatteningEnergies2: No lambda schedule specified"
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

        self.flattening_energies2_bysim = []
        self.lambda_bysim = []
        self.flattening_energies2_bylambda = {}
        for i in range(0,nlambdas):
            impact_output_file =  self.jobname + '_remd' '_%d' % i + '.out'
            datai = self.getImpactData(impact_output_file)
            fe2 = []
            lv = []
            ndata = len(datai)
            nf = len(datai[0])
            for i in range(nreject,ndata):
                lmb = datai[i][nf-4]
                eflat2 = datai[i][nf-1]
                fe2.append(eflat2)
                lv.append(lmb)
                if not lmb in self.flattening_energies2_bylambda.keys():
                    self.flattening_energies2_bylambda[lmb] = []
                self.flattening_energies2_bylambda[lmb].append(eflat2)
            self.flattening_energies2_bysim.append(fe2)
            self.lambda_bysim.append(lv)

# Saves fe/lambda trajectories for each simulation in a/lfe_trj/lfe_<sim>.dat
#
    def lfe_trj(self):
        if self.flattening_energies_bysim is None:
            self.getFlatteningEnergies()

        os.system("mkdir -p a/lfe_trj")
        for k in range(len(self.flattening_energies_bysim)):
            file = "a/lfe_trj/lfe_%d.dat" % k;
            f = open(file ,"w")
            for i in range(len(self.flattening_energies_bysim[k])):
                f.write("%f %f\n" % (self.lambda_bysim[k][i],self.flattening_energies_bysim[k][i]))
            f.close

#
# Saves fe/lambda series for each state in a/lfe_lmb/lfe_<lambda>.dat
#
    def lfe_lmb(self):
        if self.flattening_energies_bysim is None:
            self.getFlatteningEnergies()

        os.system("mkdir -p a/lfe_lmb")
        for lmb in self.binding_energies_bylambda.keys():
            file = "a/lfe_lmb/lfe_%f.dat" % lmb;
            f = open(file ,"w")
            for i in range(len(self.flattening_energies_bysim[0])):
                for k in range(len(self.flattening_energies_bysim)):
                    if self.lambda_bysim[k][i] == lmb :
                        f.write("%f %f\n" % (self.lambda_bysim[k][i],self.flattening_energies_bysim[k][i]))
            f.close

#
# Saves fe_wolam/lambda trajectories for each simulation in a/lfe_trj/lfe_wolam_<sim>.dat
#
    def lfe_wolam_trj(self):
        if self.flattening_energies2_bysim is None:
            self.getFlatteningEnergies2()

        os.system("mkdir -p a/lfe_wolam_trj")
        for k in range(len(self.flattening_energies2_bysim)):
            file = "a/lfe_wolam_trj/lfe_wolam_%d.dat" % k;
            f = open(file ,"w")
            for i in range(len(self.flattening_energies2_bysim[k])):
                f.write("%f %f\n" % (self.lambda_bysim[k][i],self.flattening_energies2_bysim[k][i]))
            f.close

#
# Saves fe/lambda series for each state in a/lfe_lmb/lfe_wolam_<lambda>.dat
#
    def lfe_wolam_lmb(self):
        if self.flattening_energies2_bysim is None:
            self.getFlatteningEnergies2()

        os.system("mkdir -p a/lfe_wolam_lmb")
        for lmb in self.binding_energies_bylambda.keys():
            file = "a/lfe_wolam_lmb/lfe_wolam_%f.dat" % lmb;
            f = open(file ,"w")
            for i in range(len(self.flattening_energies2_bysim[0])):
                for k in range(len(self.flattening_energies2_bysim)):
                    if self.lambda_bysim[k][i] == lmb :
                        f.write("%f %f\n" % (self.lambda_bysim[k][i],self.flattening_energies2_bysim[k][i]))
            f.close


# Runs MBAR using binding energies and flattening energies
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
        if self.flattening_energies2_bysim is None:
            msg = "runMBAR: Internal error: flattening energies are not defined"
            self.exit(msg)

        self.mbarpath=self.keywords.get('MBAR_PATH')
        if self.mbarpath is None:
            msg = "runMBAR: MBAR_PATH is not defined"
            self.exit(msg)
        sys.path.append(self.mbarpath)
        import MBAR

        K = len(self.lambdas)
        n = len(self.binding_energies_bysim[0])
        m = len(self.flattening_energies2_bysim[0])
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
        fe_klt = zeros([K,K,m], float64)
        fe_kt = zeros([K,m], float64)
        rp_klt = zeros([K,K,n], float64)
        for k in range(K):
            for t in range(N_k[k]):
                u_kt[k,t] = self.binding_energies_bysim[k][t]
                fe_kt[k,t] = self.flattening_energies2_bysim[k][t]
        for k in range(K):
            for t in range(N_k[k]):
                for l in range(K):
                    # Reduced energy of snapshot t from simulation k 
                    # in state l
                    u_klt[k,l,t] = beta * self.lambdas[l] * u_kt[k,t]
                    if self.lambdas[l] < 0.5:
                        fe_klt[k,l,t] = beta * (2 * self.lambdas[l]) * fe_kt[k,t]
                    if self.lambdas[l] >= 0.5:
                        fe_klt[k,l,t] = beta * (2-(2 * self.lambdas[l])) * fe_kt[k,t]
                    rp_klt[k,l,t] = u_klt[k,l,t] - fe_klt[k,l,t]

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
        mbar = MBAR.MBAR(rp_klt, N_k, relative_tolerance = 1.0e-6, verbose = True )
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
                wt = N_k[:]*exp(fk[:] - rp_klt[k,:,t])
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




##################### MAIN CODE ##########################
if __name__ == '__main__':
    
    # Setup the logger
    logger = schrodinger.utils.log.get_output_logger("bedam_analyze")

    # Parse arguments:
    usage = "%prog [options] <inputfile>"
    parser = cmdline.SingleDashOptionParser(usage)
    (options, args) = parser.parse_args(sys.argv[1:])
    
    if len(args) != 1:
        parser.error("Please specify ONE input file")
    
    commandFile = args[0]

    print ""
    print "===================================="
    print "       BEDAM Job Analysis           "
    print "===================================="
    print ""
    print "SCHRODINGER: " + os.environ['SCHRODINGER']
    print "Started at: " + str(time.asctime())
    print "Input file:", commandFile
    print ""
    sys.stdout.flush()
    
    print "Reading options"
    bedam = bedam_analyze_acflat(commandFile, options)
    print "Collecting binding energies ..."
    bedam.getBindingEnergies()
    print "Collecting flattening energies ..."
    bedam.getFlatteningEnergies()
    print "Collecting flattening energies wolam..."
    bedam.getFlatteningEnergies2()
    print "Writing binding energy histograms/trajectories and flattening energy/traj..."
    bedam.be_histograms()
    bedam.lbe_trj()
    bedam.lbe_lmb()
    bedam.lfe_trj()
    bedam.lfe_lmb()
    bedam.lfe_wolam_trj()
    bedam.lfe_wolam_lmb()
    print "Running MBAR:"
    (Deltaf_ij, dDeltaf_ij) = bedam.runMBAR()
    #Assume binding free energy is G(maxlambda)-G(minlambda)
    #typically maxlambda=1 and minlambda=0
    nf = len(Deltaf_ij[0])-1
    print "DG(binding) = %f +- %f kcal/mol" % (Deltaf_ij[0][nf], dDeltaf_ij[0][nf])
    print Deltaf_ij

#EOF

