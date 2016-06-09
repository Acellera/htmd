# Suppress the unhelpful warning  that matplotlib emits
import warnings
warnings.filterwarnings("ignore")

import matplotlib as mpl

mpl.use('Agg')

import os
import shutil
import re
import sys
from htmd.parameterization.configuration import ParameterizationConfig
from htmd.parameterization.command import Command

from htmd.parameterization.parameterization import Parameterization
import logging


# from htmdx.cli import *

def syntax():
    print("")
    print(" Acellera Small Molecule Parameterisation 2016 (c) \n")
    print("  Syntax: parameterization [args]:\n")
    print("    --input [inputfile]      : Specify the input file. Default value 'input'")
    print("    --resume [job]           : Specify name of a job to resume")
    print("    --command ([command])    : Show help for an input-file command")
    print("    --list                   : Show jobs")
    print("    --rename [molfile]       : Rename the atoms in specified file")
    print("    --torsions [molfile]     : List the soft torsions")
    print("    --delete [job]           : Delete a named job")
    print("    --license                : Show the software license")
    print("    --verbose                : Verbose logging")
    print("    --cache                  : Cache directory [~/.htmd/gaamp/]")
    print("    --help")
    print("")
    print(" No re-distribution in whole or part")
    print("")


def main_parameterize():
    import matplotlib as mpl
    mpl.use('Agg')

    input_file = None
    device = None
    verbose = False
    platform = None
    debug = False
    list_jobs = False
    rename_file = None
    job_to_delete = False
    resume = None
    cache_prefix = None
    # check_registration( product="parameterization" )

    a = 1
    try:
        while (a < len(sys.argv)):
            if (sys.argv[a] == "--rename"):
                if a < len(sys.argv) - 1:
                    rename_file = sys.argv[a + 1]
                    a = a + 1
                else:
                    raise NameError()
            elif (sys.argv[a] == "--cache"):
                if a < len(sys.argv) - 1:
                    cache_prefix = sys.argv[a + 1]
                    Parameterization.set_prefix(cache_prefix)
                    a = a + 1
                else:
                    raise NameError()

            elif (sys.argv[a] == "--torsions"):
                if a < len(sys.argv) - 1:
                    try:
                        tt = (Parameterization.listDihedrals(sys.argv[a + 1]))
                        tt = tt[1]
                        for x in tt:
                            for y in x:
                                print("%5s" % (y), end="")
                            print("")
                    except NameError as e:
                        print(str(e))
                    sys.exit(0)
                else:
                    raise NameError()

            elif (sys.argv[a] == "--input"):
                if a < len(sys.argv) - 1:
                    input_file = sys.argv[a + 1]
                    a = a + 1
                else:
                    raise NameError()
            elif (sys.argv[a] == "--resume"):
                if a < len(sys.argv) - 1:
                    resume = sys.argv[a + 1]
                    a = a + 1
                else:
                    raise NameError()
            elif sys.argv[a] == "--license":
                #				show_license( product="parameterization" )
                sys.exit(0)
            elif (sys.argv[a] == "--debug"):
                debug = True
            elif (sys.argv[a] == "--list"):
                list_jobs = True
            elif (sys.argv[a] == "--delete") and (a < len(sys.argv) - 1):
                job_to_delete = sys.argv[a + 1]
                a = a + 1
            elif (sys.argv[a] == "--verbose"):
                verbose = True
            elif (sys.argv[a] == "--command"):
                if a < len(sys.argv) - 1:
                    Command.help(sys.argv[a + 1])
                    sys.exit(1)
                    a += 1
                else:
                    Command.help(None)
                    sys.exit(1)

            elif sys.argv[a] == "--help" or sys.argv[a] == "-h":
                syntax()
                sys.exit(0)
            else:
                print("Unparsed : " + sys.argv[a])

            a += 1
    except NameError as e:
        if debug:
            raise e
        syntax()
        sys.exit(0)

    if rename_file:
        try:
            Parameterization.renameStructure(rename_file)
            print("Atoms renamed and file over-written")
            sys.exit(0)
        except (NameError, ValueError) as e:
            print("Renaming of atoms failed " + str(e))
            sys.exit(1)

    if input_file and resume:
        print("--input and --resume are mutually exclusive")
        sys.exit(1)

    if list_jobs:
        Parameterization.listJobs()
        sys.exit(0)
    if job_to_delete:
        try:
            Parameterization.deleteJob(job_to_delete)
            sys.exit(0)
        except (ValueError, NameError) as e:
            print("Unable to remove job: " + str(e))
            sys.exit(1)

    if not resume and (not input_file):
        input_file = "input"

    if input_file and not os.path.isfile(input_file):
        print("Input file not found")
        syntax()
        sys.exit(0)

    print("\n === Parameterise 2016 ===\n")
    print("      (c) Acellera ")
    try:
        if input_file:
            config = ParameterizationConfig(config=input_file, check=True)
            print(config)
        else:
            config = None
    except NameError as e:
        print("\n Failed to parse input file : \n\n  " + str(e) + "\n\n")
        if debug:
            raise e
        sys.exit(2)

    if verbose:
        logging.getLogger("htmd.parameterization").setLevel(logging.INFO)
    else:
        logging.getLogger("htmd.parameterization").setLevel(logging.WARNING)

    if debug:
        logging.getLogger("htmd.parameterization").setLevel(logging.DEBUG)
        if config:
            config.Debug = True

    try:
        p = Parameterization(config=config, name=resume)
        (fp, names) = p.plotDihedrals(show=False)
        pp = p.getParameters()
        for i in range(len(fp)):
            shutil.copyfile(fp[i], "torsion-potential-" + re.sub(" ", "-", names[i]) + ".svg")
        shutil.copyfile(pp['RTF'], "mol.rtf")
        shutil.copyfile(pp['PRM'], "mol.prm")
        shutil.copyfile(pp['PDB'], "mol.pdb")

    except (NameError, ValueError) as e:
        print("\n Failed to run parameterisation : " + str(e) + "\n")
        if debug:
            raise e
        sys.exit(3)

    except KeyboardInterrupt as e:
        print("\n Terminating at user request before parameterisation started\n")
        if debug:
            raise e
        sys.exit(5)

    sys.exit(0)


if __name__ == "__main__":
    main_parameterize()
