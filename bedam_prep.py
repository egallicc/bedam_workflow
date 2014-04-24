# python script
__doc__="""
$Revision: 0.1 $

BEDAM job preparation script

"""
# Contributors: Emilio Gallicchio, Junchao Xia

import os, sys, time
from schrodinger import structure, structureutil
from schrodinger.application import inputconfig
from schrodinger.utils import cmdline
import schrodinger.utils.log

from bedam import bedam_job

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
    bedam = bedam_job(commandFile, options)
    if bedam.impact_package is None:
        msg = "bedam_prep: No impact package specified in the input file"
        bedam.exit(msg)
   
    if  bedam.impact_package == "commercial" :
        print "Analyzing structure files ..."
        bedam.getImpactMaeFiles()
        print "Writing BEDAM restraint file ..."
        bedam.writeRestraintFile()
        print "Adding atomic restraints to receptor file ..."
        bedam.writeRecStructureFile()
    elif bedam.impact_package == "academic" :   
        print "Analyzing structure files ..."
        bedam.getDesmondDMSFiles()
        print "Writing BEDAM restraint file ..."
        bedam.writeRestraintFileFromDMS()
        print "Adding atomic restraints to receptor file ..."
        bedam.writeRecStructureDMSFile()
    elif bedam.impact_package == "wcg" :
        print "Analyzing structure files ..."
        bedam.getDesmondDMSFiles()
        print "Writing BEDAM restraint file ..."
        bedam.writeRestraintFileFromDMS()
        print "Adding atomic restraints to receptor file ..."
        bedam.writeRecStructureDMSFile()
    else: 
        msg = "bedam_prep: The impact package specified in the input file is not right"
        bedam.exit(msg)
       

    print "Writing Minimization/Thermalization input file ..."
    bedam.writeThermInputFile()
    print "Writing Remd Production input file ..."
    bedam.writeRemdInputFile()
    print "Writing Trajectory conversion input file ..."
    bedam.writeReadTrajInputFile()

    print
    print "Job preparation complete"
    print ""
    if  bedam.impact_package == "commercial" :	
    	print "To run the minimization/thermalization calculation do:"
    	print "$SCHRODINGER/impact -i <jobname>_mintherm.inp -HOST ... -LOCAL"
    	print ""
   	print "When completed run the production calculation with:"
    	print "$SCHRODINGER/impact -i <jobname>_remd.inp -HOST ... -LOCAL"
    	print ""
    	print "When completed run the BEDAM analysis with:"
    	print "$SCHRODINGER/run ~/utils/bedam_scripts/bedam_analyze.py <jobname>.inp"
    	print "  The computed binding free energy will be shown"
    	print "  Also, binding energy histograms at each lambda will be found in a/be_hist/"
    	print "  and lambda/binding energy trajectories for each replica will be found in a/lbe_trj/"
    	print ""
    	print "Optionally create conformational trajectories for each lambda with:"
    	print "$SCHRODINGER/impact -i <jobname>_readtraj.inp -HOST ... -LOCAL"
    	print " These will be named <jobname>_remd_trj_<lambda>.maegz"

#EOF
