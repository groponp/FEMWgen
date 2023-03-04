#!----------------------------------------------------------------------------------#!
#! FEMWgen  : A computational tool to generate multiple windows                     #!
#!            for free-energy calculations.                                         #!
#!                                                                                  #!
#! @authors : Ropón-Palacios G., Alexandre Suman de Araujo e Luiz Fernando Zonetti. #!
#! date     : Tue 8 Nov., 2022.                                                     #!
#! E-mail   : <asaraujo@ibilce.unesp.br>, <groponp@gmail.com>.                      #!
#! version  : 1.0b.                                                                 #!
#!----------------------------------------------------------------------------------#!


#! Call class and methods.
#!--------------------------------

import optparse
from colorama import Fore
disclaimer = """<FEMWgen v1.0b: A computational tool to generate multiple windows for free-energy calculations.>"""
parser = optparse.OptionParser(description=disclaimer)

	#!INPUTS
parser.add_option("--psff", help="PSF file name.", type=str)
parser.add_option("--pdbf", help="PDB file name.", type=str)
parser.add_option("--bin", help="Spacing for each windows.", default=2, type=int)
parser.add_option("--lresn", help="Resname lipid component into membrane. i.e. \"POPC POPG\".", type=str, action='store')
parser.add_option("--pref", help="Reference atoms from protein. i.e. \"resname TRP and name CH2\".", type=str, action='store')
parser.add_option("--lref", help="Reference atoms from lipids. \"resname POPC POPG and name P\".")
parser.add_option("--lrst", help="Serial key of atoms to restraint during equilibration. i.e. \"serial 38224 40079 39148 39811 20055 14892 19787 18877\".", type=str, action='store')
parser.add_option("--namd", help="PATH to namd executable. i.e. \"$HOME/Documents/packages/NAMD_2.14_MacOSX-x86_64-multicore/namd2\".", type=str)
parser.add_option("--nproc", help="Number of CPUs that will be use into energy minimization.", type=int)
parser.add_option("--charm", help="PATH to charmrun executable, only if you have it compile with namd. Default \"False\".", default="False", type=str)
parser.add_option("--cfw", help="Center of first windows. Default -38", default=-38, type=int)
parser.add_option("--seltxt", help="Select protein. i.e. \"protein\".", type=str)
parser.add_option("--nw", help="Number of windows. Default 60.", default=60, type=int)


	#! Only RESTART
parser.add_option("--restart", help="Only use it option if you like restart your system preparation, from checkpoint. Default \"False\".", default="False", type=str)
parser.add_option("--chkpoint", help="Pass checkpoint file, if you like start system preparation", type=str)
parser.add_option("--jsonfile", help="JSON file contaning all paramters inputs, for restart system preparation.", default="FEMWgen.json", type=str)

	#! Usage
parser.add_option("--usage", help="Use it option to print into screen, example of command line.", action="store_true" default=True, dest="usage")


opts, args = parser.parse_args()

	#! Calling procedure.
	#!-----------------------------------



if opts.usage:
	print("Common usage NAMD only:")
	print("-----------------------")
	print("python FEMWgen.py --psff=sistema.psf --pdbf=sistema.pdb --bin=2 --lresn=\"POPC POPG\" --pref=\"resname TRP and name CH2\"")
	print("--lref=\"resname POPC POPG and name P\" --lrst=\"serial 38224 40079 39148 39811 20055 14892 19787 18877\"")
	print("--namd=\"$HOME/Documents/packages/NAMD_2.14_MacOSX-x86_64-multicore/namd2\" --nproc=2 --cfw=-38 --seltxt=\"protein\" --nw=60\n")

	print("Common usage NAMD w/ CHARMM++:")
	print("------------------------------")
	print("python FEMWgen.py --psff=sistema.psf --pdbf=sistema.pdb --bin=2 --lresn=\"POPC POPG\" --pref=\"resname TRP and name CH2\"")
	print("--lref=\"resname POPC POPG and name P\" --lrst=\"serial 38224 40079 39148 39811 20055 14892 19787 18877\"")
	print("--namd=\"$HOME/Documents/packages/NAMD_2.14_MacOSX-x86_64-multicore/namd2\" --nproc=2 --charm=\"$HOME/Documents/packages/charmrun/charmrun\" --cfw=-38 --seltxt=\"protein\" --nw=60\n")

	print("Restart usage:")
	print("--------------")
	print("python FEMWgen.py --restart=\"True\" --chkpoint=FEMWgen.chkpoint --jsonfile=FEMWgen.json --cfw=-38 --seltxt=\"protein\" ")



else: 

	if opts.restart == "False":
		import FEMWgenLIB

		mkw = FEMWgenLIB.MakeWindows(psff=opts.psff, pdbf=opts.pdbf, Bin=opts.bin, lresn=opts.lresn, refprot=opts.pref, refmemb=opts.lref, rstmemb=opts.lrst,
                   namdPATH=opts.namd, nproc=opts.nproc, restart=opts.restart, charmrunPATH=opts.charm)

		FEMWgenLIB.IO.writrJsonParms(ofile="FEMWgen.json", dic=mkw.__dict__)
		mkw.iterateWindows(cfw=-38, seltxt=opts.seltxt, N=opts.nw)


	elif opts.restart == "True":
		import FEMWgenLIB
		restart = FEMWgenLIB.RestartWindows()
		windowsOK = restart.get_NotiterateWindows(opts.chkpoint)
		restart.restart(windowsOK, cfw=opts.cfw, seltxt=opts.seltxt, jsonfile=opts.jsonfile)

	else:
		print(Fore.RED +  "Error no command pass." + Fore.RESET)

