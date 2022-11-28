'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Jean Mainguy, 28 nov. 2022 
License     : MIT 
Maintainer  : jean.mainguy@inrae.fr 
Portability : POSIX


'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def init_logging(args):
    '''Initialise logging.'''

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt="%Y-%m-%dT%H:%M:%S%z")

    else:
        logging.basicConfig(format='%(asctime)s %(levelname)s - %(message)s',
                    datefmt="%Y-%m-%dT%H:%M:%S%z")

    logging.info('Program started')
    logging.info(f'command line: {" ".join(sys.argv)}', )



def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(description="...",
                            formatter_class=ArgumentDefaultsHelpFormatter)
                        

    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")


    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)

    args = parser.parse_args()
    return args


def main():
    "Orchestrate the execution of the program"

    args = parse_arguments()
    init_logging(args)


    ### check input validity 

    ### collect info on all bins sets  
        # bin set is made of multiple bins
            # bin is a list of contigs
            # Class Bin 
            ## attributes 
            #  contigs = list of contigs
            #  origin = name of the orginated bin set 
            #  id = id of the bin 
            
    ### apply binning_refiner on all bins sets pairs 
        ### parallelize by number of possible cpus ? 

    ### run checkm2 on all bins sets
        #### initialise checkm2 
            #### run prodigal or use faa input files 
                ##### input faa need faa files and gff file to have contig2faa_names
            #### run diamond on all faa 

        #### run checkm2 on all bin sets in a parallel fashion
            #### organise checkm2 temp files 
            #### run checkm2 with --resume flag

    ### consolidate bins

    ### run checkm2 on final bin sets

    ### write final bin sets 

    ### plot bins


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
