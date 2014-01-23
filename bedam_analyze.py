# python script
__doc__="""
$Revision: 0.1 $

BEDAM analysis script

"""
# Contributors: Emilio Gallicchio

import os, sys, time
from operator import itemgetter, attrgetter
from schrodinger import structure, structureutil
from schrodinger.application import inputconfig
from schrodinger.utils import cmdline
import schrodinger.utils.log

from bedam import bedam_job

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
    bedam = bedam_job(commandFile, options)
    print "Collecting binding energies ..."
    bedam.getBindingEnergies()
    print "Writing binding energy histograms/trajectories ..."
    bedam.be_histograms()
    bedam.lbe_trj()
    print "Running MBAR:"
    (Deltaf_ij, dDeltaf_ij) = bedam.runMBAR()
    #Assume binding free energy is G(maxlambda)-G(minlambda)
    #typically maxlambda=1 and minlambda=0
    nf = len(Deltaf_ij[0])-1
    print "DG(binding) = %f +- %f kcal/mol" % (Deltaf_ij[0][nf], dDeltaf_ij[0][nf])
    print Deltaf_ij

#EOF
