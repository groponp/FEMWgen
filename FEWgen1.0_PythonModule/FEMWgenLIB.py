#!----------------------------------------------------------------------------------#!
#! FEMWgen  : A computational tool to generate multiple windows                     #!
#!            for free-energy calculations.                                         #!
#!                                                                                  #!
#! @authors : Ropón-Palacios G., Alexandre Suman de Araujo e Luiz Fernando Zonetti. #!
#! date     : Tue 8 Nov., 2022.                                                     #!
#! E-mail   : <asaraujo@ibilce.unesp.br>, <groponp@gmail.com>.                      #!
#! version  : 1.0b.                                                                 #!
#!----------------------------------------------------------------------------------#!


#!---------------------------------#!
#! ChangeLogs:                     #!
#!---------------------------------#!
#! 1) This is an reimplementation from original code writed in TCL. 14:03 pm, Tue 8, 2022. 
#! 2) Add to method writeFiles, `package require psfgen` and `resetpsf`, to fix atoms incongruency. 19:30 pm, Wed 9 Nov, 2022.
#! 3) Add class Restart. 20:20, Mon 21, Nov, 2022.



#! Import Libraries
#!----------------------
from vmd import evaltcl 
from vmd import molecule
import os 
from colorama import Fore 
import time
import glob
import json

#! Create Class
#!-----------------------

class IO:

	def __init__(self):
		pass

	def loadFiles(self, psff, pdbf):
		molid = molecule.load("psf", psff)
		evaltcl("mol addfile  {} wait for all {}".format(pdbf, molid))
		print(Fore.RED + "[INFO   ] Loading PSF-{} and PDB-{} to molid {}".format(psff, pdbf, molid) + Fore.RESET)
		return molid

	def writeFiles(self, id, ofile, selxtAtoms, mode="both"):

		if mode == "both":
			sel = evaltcl("atomselect {} \"{}\"".format(id, selxtAtoms))
			evaltcl("package require psfgen")
			evaltcl("resetpsf")
			evaltcl("{} writepsf ./{}.psf".format(sel, ofile))
			evaltcl("{} writepdb ./{}.pdb".format(sel, ofile))
			print(Fore.RED + "[INFO   ] Following files [{}.psf, {}.pdb] was writed.".format(ofile, ofile) + Fore.RESET)

		elif mode == "pdb":
			sel = evaltcl("atomselect {} \"{}\"".format(id, selxtAtoms))
			evaltcl("{} writepdb ./{}.pdb".format(sel, ofile))
			print(Fore.RED + "[INFO   ] Following file {}.pdb was writed.".format(ofile) + Fore.RESET)

		elif mode == "psf":
			sel = evaltcl("atomselect {} \"{}\"".format(id, selxtAtoms))
			evaltcl("{} writepsf ./{}.psf".format(sel, ofile))
			print(Fore.RED + "[INFO   ] Following file {}.psf was writed.".format(ofile) + Fore.RESET)

		else:
			print(Fore.RED + "[ERROR.  ] You haven't select  mode." + Fore.RESET)


	def makeFolders(self, N=60, ofolder="w"):
		c = 1 
		while c <= N:

			if os.path.exists("{}_{}".format(ofolder, c)):
				#! All this file will be int configs folder: em_janela.conf em.conf md_eq.conf md_abf.conf md_abf.in par_all27_prot_na.prm par_all36_lipid.prm. 
				os.system("cp configs/* {}_{}/".format(ofolder, c)) 
				print(Fore.RED + "[INFO   ] Folder {}_{} already created, only copy files.".format(ofolder, c) + Fore.RESET)
				c += 1

			else: 
				os.mkdir("{}_{}".format(ofolder, c))
				os.system("cp configs/* {}_{}/".format(ofolder, c)) 
				print(Fore.RED + "[INFO   ] Making folder {}_{} and copy files.".format(ofolder, c) + Fore.RESET)
				c += 1

	def selector(self, id, seltxt):
		sel = evaltcl("atomselect {} \"{}\"".format(id, seltxt))
		print(Fore.RED + "[INFO   ] Selecting {} Atoms.".format(seltxt) + Fore.RESET)
		return sel 

	def checkTOP(self):
		evaltcl("mol delete top")           #! Check out to not keeping any mol loaded. 
		list  = evaltcl("mol list")
		check = list.split()[0]
		check = check.replace("ERROR)", "True")
		return check 

	def writrJsonParms(ofile, dic):
		dictionary = dic 
		json_object = json.dumps(dictionary, indent=10)
		with open (ofile, "w") as outfile:
			outfile.write(json_object)



class GeomTrans: 

	def __init__(self):
		pass

	def calc_dz(self, cfw, ref_memb, ref_prot):
		dz = float(evaltcl("expr ({} - [lindex [measure center {}] 2]) - ([lindex [measure center {}] 2] - [lindex [measure center {}] 2])".format(cfw,ref_memb, ref_prot, ref_memb)))
		return dz

	def translateAtoms(self, dz, mobilAtoms):
		vec = evaltcl("lappend moveprot 0 0 {}".format(dz))
		evaltcl("{} moveby [list {}]".format(mobilAtoms, vec))
		evaltcl("unset moveprot")
		print(Fore.RED + "[INFO   ] Translate {} Atoms to the vector {}.".format(mobilAtoms, vec) + Fore.RESET)

	def translatePeptide(self, id, dz, cfw, Bin, c):
		cf = evaltcl("expr {}+({}-1)*{}".format(cfw, c, self.Bin))
		dz = evaltcl("expr {} - {}".format(cf, cfw))
		prot = evaltcl("atomselect {} protein".format(id))
		moveprot = evaltcl("lappend moveprot 0 0 {}".format(dz))
		evaltcl("{} moveby [list {}]".format(prot, moveprot))
		evaltcl("unset moveprot")
		print(Fore.RED + "[INFO   ] Translate peptide to the vector 0 0 {}.".format(dz) + Fore.RESET)
		return int(cf) 

	def zminmax(self, seltxt):
		min = float(evaltcl("lindex [lindex [measure minmax {}] 0] 2".format(seltxt)))
		max = float(evaltcl("lindex [lindex [measure minmax {}] 1] 2".format(seltxt)))
		return min, max  

	def meanZ_position(self, id, lresn):
		pmin = evaltcl("atomselect {} \"resname {} and name P and z < 0\"".format(id, lresn))
		pmax = evaltcl("atomselect {} \" resname {} and name P and z > 0\"".format(id, lresn))
		pmin_z = float(evaltcl("lindex [measure center {}] 2".format(pmin)))
		pmax_z = float(evaltcl("lindex [measure center {}] 2".format(pmax)))
		return pmin_z, pmax_z

	def distanceZ(self, pmin_z, pmax_z, min, max):
		z1 = float(evaltcl("expr {}-{}".format(pmin_z,max)))
		z2 = float(evaltcl("expr {}-{}".format(pmax_z,min)))
		return z1,z2

	def rmsd(self, pos_x_ref, pos_x_tmp):
		j = 0
		rmsd = 0
		while j < int(evaltcl("llength {}".format(pos_x_ref))):
			rmsd = float(evaltcl("set rmsd [expr {} + (([lindex {} {}]-[lindex {} {}])*([lindex {} {}]-[lindex {} {}]))]".format(rmsd, pos_x_tmp, j, pos_x_ref, 
				j, pos_x_tmp, j, pos_x_ref, j)))
			j += 1

		rmsd = float(evaltcl("set rmsd [expr sqrt($rmsd/[llength {}])]".format(pos_x_ref))) 
		return rmsd


class MakeWindows:

	def __init__(self,psff, pdbf, Bin, lresn, refprot, refmemb, rstmemb, namdPATH, nproc, restart = "False", charmrunPATH="False"):
		"""
		INPUTS:
		-----------
		A) <psff(str)>       : This' the name of the PSF file that contains the topology. 
		B) <pdbf(str)>       : This' the name of the PDB file that contains the coordinates. 
		C) <Bin(int)>        : Spacing of winwonds. Example: 2.
		D) <lresn(str)>      : resname of lipid components into membrane. Example 'resname POPC POPG'.
		E) <refprot(str)>    : Reference protein selection.
		F) <refmemb(str)>    : Reference membrane selection.
		G) <rstmemb(str)>    : Restraint selection for membrane equilibration. 
		H) <charmrunPATH(str)> : PATH to executable of charmmrun. If charmmrun is into PATH only pass: charmrun
		I) <namdPATH(str)>   : PATH to excutable of namd.
		J) <nprocs(int)>     : Number of CPus to be used with NAMD.
		#**K) <gpuList(str)>   : gpuLisT is ide from gpu. Example in string format "0,1,2". 

		OUTPUTS:
		------------
		A) Return each argument as an attribute, to use in others methods.

		"""
		self.psff      = str(psff)
		self.pdbf      = str(pdbf)
		#self.seltxt    = str(seltxt)
		self.Bin       = int(Bin)
		self.lresn     = str(lresn)
		self.refprot   = str(refprot)
		self.refmemb   = str(refmemb)
		self.rstmemb   = str(rstmemb)
		self.charmrunPATH = str(charmrunPATH)
		self.namdPATH  = str(namdPATH)
		self.nproc     = int(nproc)
		self.restart   = str(restart)
		#self.gpuList   = gpuList



	if os.path.exists("FEMWgen.chkpoint"): 
		print(Fore.RED +  "[INFO   ] Please use your check point file to continue your system preparation." + Fore.RESET) 

	else: 

		chkpoint = open("FEMWgen.chkpoint", "w")
		chkpoint.write("#! It check point was write by FEMWgen v1.0b.\n")
		chkpoint.write("#! @uthors: Ropón-Palacios G., Alexandre Suman de Araujo e Luiz Fernando Zonetti.\n")
		chkpoint.write("#! contact to: <asaraujo@ibilce.unesp.br>, <groponp@gmail.com>.\n\n")
		chkpoint.write("-----------------------------------------------------------------------------------\n")
		chkpoint.write("NAME\t\tSTATUS\t\tTIME\n")
		chkpoint.write("-----------------------------------------------------------------------------------\n")
		chkpoint.close()
		
	def iterateWindows(self, cfw, seltxt, N=60, restartfrom="None"):
		"""
		INPUTS:
		-----------
		A) <cfw(int)>    : Center for first windows.
		B) <N(int)>      : Number of windows that user like to create. Default is 60.
		C) <seltxt(str)> : A text selecting the atoms to be moved in VMD syntax. 

		OUTPUTS:
		------------

		"""

		if self.restart == "False":
			IO.makeFolders(self, N=N, ofolder="w")

		if restartfrom == "None":
			c = 1

		else: 
			c = int(restartfrom.split("_")[1])

		while c <= N:

			check = IO.checkTOP(self)

			if check == "True":


				molid1   = IO.loadFiles(self, self.psff, self.pdbf) 
				prot     = IO.selector(self, molid1, seltxt)
				ref_memb = IO.selector(self, molid1, self.refmemb)
				ref_prot = IO.selector(self, molid1, self.refprot)

				
				#! Calculate distance in z (dz).
				#!----------------------------------------------------------
				dz = GeomTrans.calc_dz(self, cfw, ref_memb, ref_prot)
				GeomTrans.translateAtoms(self, dz, mobilAtoms=prot)

				#! Create Folders and copy files.
				#!----------------------------------------------------------
				#!print(Fore.RED + "[INFO   ] Entry to windows: {}.".format(c) + Fore.RESET)
				#!if os.path.isfile("w_{}/sistema.psf".format(c)) and os.path.isfile("w_{}/sistema.pdb".format(c)) : os.system("rm w_{}/sistema.p*".format(c))
				
				if restartfrom == "None":
					print(Fore.RED + "[INFO   ] Entry to windows: {}.".format(c) + Fore.RESET)
					if os.path.isfile("w_{}/sistema.psf".format(c)) and os.path.isfile("w_{}/sistema.pdb".format(c)) : os.system("rm w_{}/sistema.p*".format(c))
					os.chdir("w_{}".format(c))
					#os.system("cd w_{}".format(c))

				else: 
					print(Fore.RED + "[INFO   ] Entry to windows: {}.".format(restartfrom.split("_")[1]) + Fore.RESET)
					x = restartfrom.split("_")[1]
					#os.system("cd w_{}".format(x))
					os.chdir("w_{}".format(x))
					if os.path.isfile("w_{}/sistema.psf".format(x)) and os.path.isfile("w_{}/sistema.pdb".format(x)) : os.system("rm w_{}/sistema.p*".format(x))
					

				time.sleep(5)
				print(Fore.RED + "[INFO   ] Current dir: {}".format(os.getcwd()) + Fore.RESET)


				#!  Translate peptide. 
				#!---------------------------------------------------------
				cf = GeomTrans.translatePeptide(self, molid1, dz, cfw, self.Bin, c)

				#! Write file PSF and PDB after translate.
				#!---------------------------------------------------------
				IO.writeFiles(self, molid1, ofile=self.pdbf.split(".")[0], selxtAtoms="not same residue as water within 2 of protein", mode="both") 

				#! load MOLs file again and calculate zmin and zmax.
				#!--------------------------------------------------------
				molid2 = IO.loadFiles(self, self.psff, self.pdbf)
				prot2  = IO.selector(self, molid2, "protein")      
				zmin,zmax = GeomTrans.zminmax(self, prot2)  

				#! Calculate z-position, distance.
				#!-------------------------------------------------------
				pminz, pmaxz = GeomTrans.meanZ_position(self, molid2, self.lresn)
				z1,z2 = GeomTrans.distanceZ(self, pminz, pmaxz, zmin, zmax)
				#print(z1,z2)
				#time.sleep(60)


				#! If protein're into inside of membrane make open hole, if not make 
				#! restraint and update .in to next windows.
				#!------------------------------------------------------
				if z1 <= 2 and z2 >= -2:
					#print("{} <= 2  and {} >= -2".format(z1, z2))
					#time.sleep(60)

					#! Build selection the atoms of reference and get x coordinate value.
					#!---------------------------------------------------------------------
					print(Fore.RED + "[INFO   ] Peptide entry into hole of membrane.")
					atom_ref = IO.selector(self, molid2, self.rstmemb)
					pos_x_ref = evaltcl("list [{} get x]".format(atom_ref))

					#! Centrer and Lenght the protein in X axis.
					#!-----------------------------------------------------------------------
					center_x_prot = float(evaltcl("lindex [measure center {} ] 0".format(prot2)))
					lenght_x_prot = float(evaltcl("expr [lindex [lindex [measure minmax {}] 1] 0] - [lindex [lindex [measure minmax {}] 0] 0]".format(prot2, prot2))) 
					print(Fore.RED + "[INFO   ] Defining the bands." + Fore.RESET)

					#! Defining the bands to open hole into membrane.
					#------------------------------------------------------------------------
					memb = IO.selector(self, molid2, "resname "+self.lresn)
					residue_mol = int(evaltcl("list [lindex [{} get residue] 0]".format(memb)))

					#! Iterate each element da selection the membrane and check out to from that face is.
					#!------------------------------------------------------------------------
					rmol = IO.selector(self, molid2, "residue "+str(residue_mol)) 
					control = float(evaltcl("lindex [measure center {}] 0".format(rmol))) 
					
					if control < center_x_prot:
						tmp = evaltcl("lappend residue_banda_inf {}".format(residue_mol))      #! below.
						residue_banda_inf = tmp

					else: 
						tmp = evaltcl("lappend residue_banda_sup {}".format(residue_mol))   #! Up
						residue_banda_sup = tmp

					membList = evaltcl("{} get residue".format(memb))
					membList = membList.split()

					#! set residue_banda_inf to x, to avoid `error : local variable referenced before assgiment`
					residue_banda_inf = "x"

					#! Storing residues to upper o below membrane face.
					#!-------------------------------------------------------------------------
					for itemp in membList:
						if itemp != residue_mol:
							residue_mol = itemp 
							control = float(evaltcl("lindex [measure center {}] 0".format(rmol)))
							if control < center_x_prot and  control != "x":
								tmp = evaltcl("lappend residue_banda_inf {}".format(residue_mol)) 
								residue_banda_inf = tmp
								print(Fore.RED + "[INFO   ] Append to List below residue {}".format(residue_mol) + Fore.RESET)

							else:
								tmp = evaltcl("lappend residue_banda_sup  {}".format(residue_mol))
								residue_banda_sup = tmp
								print(Fore.RED + "[INFO   ] Append to List upper residue {}".format(residue_mol) + Fore.RESET)

					#! Setting to that face the mebrane is residue_mol.
					#!------------------------------------------------------------------------
					band_upper = IO.selector(self, molid2, "residue "+str(residue_banda_sup))
					if residue_banda_inf != "x": band_below = IO.selector(self, molid2, "residue "+str(residue_banda_inf))
					band_below = "x"

					#! Performing open the hole into membrane
					#!-------------------------------------------------------------------------
					print(Fore.RED + "[INFO   ] Making hole into membrane." + Fore.RESET)

					evaltcl("lappend d_sup [expr {}/2] 0 0".format(lenght_x_prot))
					evaltcl("lappend d_inf [expr (-1)*{}/2] 0 0".format(lenght_x_prot))
					evaltcl("{} moveby $d_sup".format(band_upper))
					if band_below != "x": evaltcl("{} moveby $d_inf".format(band_below))
					evaltcl("unset d_sup")
					evaltcl("unset d_inf")

					#! Select all and escreve PDB and Fix Atoms.
					#!-------------------------------------------------------------------------
					all = IO.selector(self, molid2, "all")
					IO.writeFiles(self, molid2, ofile=self.psff.split(".")[0], selxtAtoms="all", mode="pdb")
					evaltcl("{} set beta 0".format(all)) 
					fix = IO.selector(self, molid2, "protein")
					evaltcl("{} set beta 1".format(fix))
					IO.writeFiles(self, molid2, ofile="myfixedatoms", selxtAtoms="all", mode="pdb") 

					print(Fore.RED + "[INFO   ] Start to close hole into membrane." + Fore.RESET)

					#! copy temp files to tmp folder.
					#!-------------------------------------------------------------------------
					#IO.makeFolders(self, N=1, ofolder="tmp")
					if os.path.exists("tmp"):
						print(Fore.RED + "[INFO   ] tmp folder already exists." + Fore.RESET)
					
					else: 
							os.system("mkdir tmp/")
					
					os.system("cp {} {} myfixedatoms.pdb em_janela.conf tmp".format(self.psff, self.pdbf))
					os.chdir("tmp") 

					#! Loops of choice. 
					#!------------------------------------------------------------------------
					gatilho  = 1
					rmsd_ant = 10000.0
					while gatilho == 1:
						print(Fore.RED + "[INFO   ] Entry into loop of closing." + Fore.RESET)

						charmrun_exec = self.charmrunPATH
						namd_exec = self.namdPATH
						charm = self.charmrunPATH
						nproc = self.nproc

						if charm != "False":
							print(Fore.RED + "[INFO   ] Performing Energy minimization. Check out min.log to view progress." + Fore.RESET)
							os.system("{} +p{} {} {} > min.log ".format(charmrun_exec, nproc, namd_exec, "em_janela.conf"))
							

						else:
							print(Fore.RED + "[INFO   ] Performing Energy minimization. Check out min.log to view progress." + Fore.RESET)
							os.system("{} +p{} {} > min.log ".format(namd_exec, nproc, "em_janela.conf"))
							

						#! Load em binary coord
						#!-------------------------------------------------------------------
						molid3 = IO.loadFiles(self, self.psff, "em.coor")
						all3   = IO.selector(self, molid3, "all")
						atom_ref3 = IO.selector(self, molid3, self.rstmemb)
						pos_x_tmp = evaltcl("list [{} get x]".format(atom_ref3))

						#! Evaluate RMSD.
						#!--------------------------------------------------------------------
						rmsd = GeomTrans.rmsd(self, pos_x_ref, pos_x_tmp)
						print(Fore.RED + "[INFO   ] Current rmsd = {} and rmsd_ant = {}".format(rmsd, rmsd_ant) + Fore.RESET)
						if rmsd > rmsd_ant:
							gatilho = 0
							os.system("cp system_ant.pdb {}".format(self.pdbf))
							os.system("echo {} >> rmsd.dat".format(rmsd))
							#! add check point.
							if restartfrom == "None":
								t = time.localtime()
								ctime = time.strftime("%m/%d/%Y, %H:%M:%S", t)
								chkpoint = open("../../FEMWgen.chkpoint", "a")
								chkpoint.write("w_{}\t\tOK\t\t{}\n".format(c, ctime))
								chkpoint.close()

							else:
								t = time.localtime()
								ctime = time.strftime("%m/%d/%Y, %H:%M:%S", t)
								chkpoint = open("../../FEMWgen.chkpoint", "a")
								chkpoint.write("w_{}\t\tOK\t\t{}\n".format(restartfrom.split("_")[1], ctime))
								chkpoint.close()


						else:
							rmsd_ant = rmsd 
							os.system("echo {} >> rmsd.dat".format(rmsd))

						#! Save the minimized system.
						#!--------------------------------------------------------------------
						IO.writeFiles(self, molid3, ofile="system_ant", selxtAtoms="all", mode="pdb")

						#! Make again bands to be moved.
						#!--------------------------------------------------------------------
						banda_sup = evaltcl("atomselect {}  \"residue {}\"".format(molid3, residue_banda_sup))
						banda_inf = "x"
						if residue_banda_inf != "x": banda_inf = evaltcl("atomselect {}  \"residue {}\"".format(molid3, residue_banda_inf))
						evaltcl("{} moveby [list {}]".format(banda_sup, "-0.2 0 0"))
						if banda_inf != "x": evaltcl("{} moveby [list {}]".format(banda_inf, "0.2 0 0"))

						#! Write PDB from membrane choice.
						#!---------------------------------------------------------------------
						all3 = IO.selector(self, molid3, "all")
						fix = IO.selector(self, molid3, "protein")

						IO.writeFiles(self, molid3, ofile=self.pdbf.split(".")[0], selxtAtoms="all", mode="pdb")
						evaltcl("{} set beta 0".format(all3))
						evaltcl("{} set beta 1".format(fix))
						IO.writeFiles(self, molid3, ofile="myfixedatoms", selxtAtoms="all", mode="pdb")

					#! Copy PDB to folder before.
					#!--------------------------------------------------------------------------
					os.system("cp ../{} ../system_orig.pdb".format(self.pdbf))
					os.system("cp {} rmsd.dat ../".format(self.pdbf))
					os.chdir("../")
					#os.system("rm -rf tmp/")


				print(Fore.RED + "[INFO   ] Pass check out membrane hole." + Fore.RESET)

				#! Make fixing Atoms.
				molid4 = IO.loadFiles(self, self.psff, self.pdbf)
				all4 = IO.selector(self, molid4, "all")
				evaltcl("{} set beta 0".format(all4))
				fix4 = IO.selector(self, molid4, "protein")
				evaltcl("{} set beta 1".format(fix4))

				IO.writeFiles(self, molid4, ofile="myfixedatoms", selxtAtoms="all", mode="pdb")

				#! Make restraint to membrane P atoms.
				all4 = IO.selector(self, molid4, "all")
				evaltcl("{} set beta 0".format(all4))
				lipidHead = IO.selector(self, molid4, self.rstmemb)
				evaltcl("{} set beta 1".format(lipidHead))
				IO.writeFiles(self, molid4, ofile="restraint", selxtAtoms="all", mode="pdb")


				#! calculate Start and End of windows.
				ij = cf - (self.Bin)/2  #! Start.
				fj = cf + (self.Bin)/2  #! End.

				 

				print(Fore.RED + "[INFO   ] ij {}.".format(ij) + Fore.RESET)
				print(Fore.RED + "[INFO   ] fj {}.".format(fj) + Fore.RESET)


				#! Paste serial of atoms to reference int ABF colvar.
				main_serial = evaltcl("[atomselect {} \"{}\"] get serial".format(molid4, self.refprot))
				ref_serial = evaltcl("[atomselect {} \"{}\"] get serial".format(molid4, self.refmemb))

				print(Fore.RED + "[INFO   ] main_serial {}.".format(main_serial) + Fore.RESET)
				print(Fore.RED + "[INFO   ] ref_serial {}.".format(ref_serial) + Fore.RESET)

				#! Modified abf.in file. 
				if os.path.isfile("md_abf.in"):
					print(Fore.RED + "[INFO   ] Modified md_abf.in file.".format(c) + Fore.RESET)
					time.sleep(5)

					os.system("cat md_abf.in | sed s/SS/'{}'/ > tmp1.in".format(fj))
					os.system("cat tmp1.in | sed s/II/'{}'/ > tmp2.in".format(ij))
					os.system("cat tmp2.in | sed s/MMM/'{}'/ > tmp3.in".format(main_serial))
					os.system("cat tmp3.in | sed s/RRR/'{}'/ > md_abf.in".format(ref_serial))

				else:
					print(Fore.RED + "[NFO   ] md_abf file not exists." + Fore.RESET)
					time.sleep(5)

				#! add checkpoint.
				if restartfrom == "None":
					t = time.localtime()
					ctime = time.strftime("%m/%d/%Y, %H:%M:%S", t)
					chkpoint = open("../FEMWgen.chkpoint", "a")
					chkpoint.write("w_{}\t\tOK\t\t{}\n".format(c, ctime))
					chkpoint.close()

				else:
					t = time.localtime()
					ctime = time.strftime("%m/%d/%Y, %H:%M:%S", t)
					chkpoint = open("../FEMWgen.chkpoint", "a")
					chkpoint.write("w_{}\t\tOK\t\t{}\n".format(restartfrom.split("_")[1], ctime))
					chkpoint.close()

				if restartfrom == "None":
					print(Fore.RED + "[INFO   ] Leaving windows: {}".format(c) + Fore.RESET)
					os.chdir("../")
					time.sleep(5)

				else: 
					print(Fore.RED + "[INFO   ] Leaving windows: {}".format(restartfrom.split("_")[1]) + Fore.RESET)
				
					os.chdir("../")
					time.sleep(5)

				if restartfrom == "None":
					c += 1

				else: 
					#c = int(restartfrom.split("_")[1])
					c += 1
			
		print(Fore.RED + "[DEBUG   ] Finishing with all windows ({}) iterates using FEMWGen v1.0b".format(N) + Fore.RESET)
		
		if restartfrom == "None":
			evaltcl("quit")
	#chkpoint.close()


class RestartWindows:

	def __init__(self):
		pass

	def get_NotiterateWindows(self, chkpoint):
		f = open(chkpoint, "r")
		lines = f.readlines()
		windowsOK = []

		for l in lines :
			#print(l)
			if  l.startswith("w") and l.split()[1] == "OK":
				windowsOK.append(l.split()[0]) #! append name of windows i.e.`w_1`. 

		return windowsOK

	def restart(self, windowsOK, cfw, seltxt, jsonfile):

		dirs = glob.glob("w_*")
		x = cfw
		notIterates = []

		for element in dirs:
			if element not in windowsOK:
				notIterates.append(element)

		print(notIterates)

		for w in notIterates:
			print(Fore.RED + "[INFO ] Restart windows from w_{}".format(w.split("_")[1]) + Fore.RESET)

			
			f = open(jsonfile)

			if w != notIterates[:-1]:
				data = json.load(f)

				mkw = MakeWindows(data["psff"], data["pdbf"], data["Bin"], data["lresn"], data["refprot"], data["refmemb"], data["rstmemb"], data["namdPATH"], data["nproc"], restart="True", charmrunPATH=data["charmrunPATH"])
				num = int(w.split("_")[1])
				mkw.iterateWindows(x, seltxt, N=num, restartfrom=w)

			else: 
				data = json.load(f)

				mkw = MakeWindows(data["psff"], data["pdbf"], data["Bin"], data["lresn"], data["refprot"], data["refmemb"], data["rstmemb"], data["namdPATH"], data["nproc"], restart="True", charmrunPATH=data["charmrunPATH"])
				num = int(w.split("_")[1])
				mkw.iterateWindows(x, seltxt, N=num, restartfrom=w)

				evaltcl("quit")
						






































