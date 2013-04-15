# Licensed under the terms of the Apache Software License 2.0:
#  http://www.apache.org/licenses/LICENSE-2.0

"""
The numpy library is required. But he pylab library is 
optional, if it is not present no images will be generated
and the information will be saved in a .txt file
"""
figures = True
import sys
try: import numpy as np
except ImportError:
	print "You do not have installed the numpy module."
	sys.exit()
try: import pylab as pl
except ImportError:
	print "You do not have installed the pylab module."
	figures = False
"""
The following libraries are installed when intalling python as part
of the Standard Library
"""
import argparse
import subprocess as subp
import glob
import math

def profasi(data):
	"""
	This function implements the default regions defined in the program PROFASI.
	"""
	chain_structure = []
	for element in data:
		if element[0] > -90 and element[0] < 30: 
			if element[1] > -77 and element[1] < -17:
				state = "alpha"
			else: state = "undefined"
		elif element[0] > -150 and element[0] < -90:
			if element[1] > 90 and element[1] < 150:
				state = "beta"
			else: state = "undefined"
		else: state = "undefined"
		chain_structure.append(state)
	return chain_structure

def custom(data):
	"""
	This function implements the region defined by the user.
	"""
	chain_structure = []
	for element in data:
		state = "undefined"
		if element[0] > -int(args.customized_region[1]) and element[0] < int(args.customized_region[2]): 
			if element[1] > int(args.customized_region[3]) and element[1] < int(args.customized_region[4]):
				state = "%s"%args.customized_region[0]
		chain_structure.append(state)
	return chain_structure

def pappu (data):
	"""
	Regions defined at the article:
	- Net charge per residue modulates conformational ensembles of
		intrinsically disordered proteins
		Albert H. Mao, Scott L. Crick, Andreas Vitalis,
		 Caitlin L. Chicoine, and Rohit V. Pappu
	  PNAS, 2010, vol. 107, no 18, 8183-8188
	  at the section: Supporting information 
	  image Secondary structure definitons.
	"""
	chain_structure = []
	for element in data:
		state = "undefined"
		if element[0] >= 170 and element[0] <= 180:
			if element[1] >= 130 and element[1]<=140: state= "beta"
			elif element[1] >= 160 and element[1]<=170: state= "beta"
		elif element[0]<= -30 and element[0]>-40:
			if element[1]>=-60 and element[1]<=-20: state = "alpha"
		elif element[0] <= -40 and element[0] > -50:
			if element[1]>=-70 and element[1]<=-10: state = "alpha"
		elif element[0] <= -50 and element[0] > -60:
			if element[1] >= -70 and element[1] <= 0: state = "alpha"
			elif element[1] >= 130 and element[1] <= 150: state = "PPII"
		elif element[0] <= -60 and element[0] > -70:
			if element[1] >= -70 and element[1] <= 10: state = "alpha"
			elif element[1] >=130 and element[1] <= 180: state = "PPII"
			elif element[1] >=-180 and element[1] <= -160: state = "PPII"
		elif element[0] <= -70 and element[0] > -80:
			if element[1] >= -60 and element[1] <= 20: state = "alpha"
			elif element[1] >=120 and element[1] <= 180: state = "PPII"
			elif element[1] >=-180 and element[1] <= -160: state = "PPII"
		elif element[0] <= -80 and element[0] > -90:
			if element[1] >= -50 and element[1] <= 20: state = "alpha"
			elif element[1] >=110 and element[1] <= 180: state = "PPII"
			elif element[1] >=-180 and element[1] <= -160: state = "PPII"
		elif element[0] <= -90 and element[0] > -100:
			if element[1] >= -40 and element[1] <= 20: state = "alpha"
			elif element[1] >=110 and element[1] <= 180: state = "PPII"
			elif element[1] >=-180 and element[1] <= -160: state = "PPII"
		elif element[0] <= -100 and element[0] > -110:
			if element[1] >=110 and element[1] <= 180: state = "PPII"
			elif element[1] >=-180 and element[1] <= -160: state = "PPII"
		elif element[0] <= -110 and element[0] > -120:
			if element[1] >=120 and element[1] <= 180: state = "PPII"
			elif element[1] >=-180 and element[1] <= -160: state = "PPII"
		elif element[0] <= -120 and element[0] > -130:
			if element[1] >=110 and element[1] <= 180: state = "beta"
			elif element[1] >=-180 and element[1] <= -160: state = "beta"
		elif element[0] <= -130 and element[0] > -140:
			if element[1] >=100 and element[1] <= 180: state = "beta"
			elif element[1] >=-180 and element[1] <= -160: state = "beta"
		elif element[0] <= -140 and element[0] > -150:
			if element[1] >=100 and element[1] <= 180: state = "beta"
			elif element[1] >=-180 and element[1] <= -160: state = "beta"
		elif element[0] <= -150 and element[0] > -160:
			if element[1] >=100 and element[1] <= 180: state = "beta"
			elif element[1] >=-180 and element[1] <= -150: state = "beta"
		elif element[0] <= -160 and element[0] > -170:
			if element[1] >=110 and element[1] <= 180: state = "beta"
			elif element[1] >=-180 and element[1] <= -150: state = "beta"
		elif element[0] <= -170 and element[0] >= -180:
			if element[1] >=120 and element[1] <= 180: state = "beta"
			elif element[1] >=-180 and element[1] <= -160: state = "beta"
		chain_structure.append(state)
	return np.asarray(chain_structure)

def blackledge(data):
	"""
	Regions defined at the article:
	- Mapping the Potential Energy Landscape of Intrinsically Disordered
	  Proteins at Amino Acid Resolution
		Valery Ozenne, Robert Schneider, Mingxi Yao, Jie-rong Huang,
		Loic Salmon, Markus Zweckstetter,Malene Ringkjobing Jensen,
		and Martin Blackledge
	  J.Am.Chem.Soc, 2012, 134, 15138-15148 
	in Figure 2 page 15140.
	"""
	chain_structure = []
	for element in data:
		if element[0] > 0:
			chain_structure.append("destrogira")
		elif element[0] <0:
			if element[1]>-120 and element[1]<50:
				chain_structure.append("alpha")
			else:
				if element[0] >-90: chain_structure.append("PPII")
				else: chain_structure.append("beta")
	return np.asarray(chain_structure)

def stride(name):
	"""
	This function uses the independent program STRIDE
	to calculate the secondary structure
	"""
	call = "stride "+name+".pdb -f"+name+".out"
	call = call.split()
	subp.call(call)
	filename = open(name+".out", "r")
	lines = [line.split() for line in filename if line != " \n"]
	chain_structure = [line[5] for line in lines if line[0] == "ASG"]
	del chain_structure[0]
	del chain_structure[-1]
	for i in range(len(chain_structure)):
		if chain_structure[i] == "H": chain_structure[i] = "alpha"
	for i in range(len(chain_structure)):
		if chain_structure[i] == "E": chain_structure[i] = "beta"
	cleaning = ("rm "+name+".out")
	cleaning = cleaning.split()
	subp.call(cleaning)
	return chain_structure

def count (data, struct):
	"""
	This function count how many aminoacids are in a given structure in a row.
	"""
	mat = np.zeros([len(data), len(data)+1])
	a=[]
	for i in range(len(data)):
		if data[i]!= struct and a:
			mat[a,len(a)] = 1.
			a = []
			found = True
		else:
			if data[i] == struct:
				a.append(i)
	if a: mat[a,len(a)] = 1.
	for resi in mat:
		if all(resi[1:]==0.): resi[0]=1.
	return mat

def images (percentages, structure):
	"""
	This function generates all the images except the one with multiple temperatures.
	"""
	fig = pl.figure(structure)
	ax = fig.add_subplot(111)
	cs = ax.matshow(percentages[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]], cmap = args.cm)
	if args.rgc: cs.set_clim(args.rgc[0],args.rgc[1])
	fig.colorbar(cs)
	(ydim, xdim) =percentages[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]].shape
	starting_residue=args.residues[0]+2
	xticks_spacing=int(ax.get_xticks()[1]-ax.get_xticks()[0])
	yticks_spacing=int(ax.get_yticks()[1]-ax.get_yticks()[0])
	xlabels=["",]+[x for x in range(1,xdim+1, xticks_spacing)]+["",]
	ylabels=["",]+[y+starting_residue for y in range(0,ydim, yticks_spacing)]+["",]
	ax.set_xticklabels(xlabels)
	ax.set_yticklabels(ylabels)
	ax.grid()
	fig.show()
	if args.save_figure:
		pl.savefig("%s-estructura-%s-%s definition"%(args.save_figure[0],structure,args.structure_definition))
	return

def degrees(rad_angle) :
	"""Converts any angle in radians to degrees.

	If the input is None, the it returns None.
	For numerical input, the output is mapped to [-180,180]
	"""
	if rad_angle is None :
		return None
	angle = rad_angle * 180 / math.pi
	while angle > 180 :
		angle = angle - 360
	while angle < -180 :
		angle = angle + 360
	return angle

def angles_calcul (pdb_code):
	import Bio.PDB
	structure = Bio.PDB.PDBParser().get_structure(pdb_code, "%s.pdb" % pdb_code)
	angles = []
	for model in structure :
		for chain in model :
			#print "Chain %s" % str(chain.id)
			polypeptides = Bio.PDB.CaPPBuilder().build_peptides(chain)
			for poly_index, poly in enumerate(polypeptides) :
				phi_psi = poly.get_phi_psi_list()
				for res_index, residue in enumerate(poly) :
					phi, psi = phi_psi[res_index]
					if phi and psi :
						#Don't write output when missing an angle
						angles.append([degrees(phi), degrees(psi)])
	return np.asarray(angles)

def pdb_npy (folder):
	filelist = glob.glob(folder + '*.pdb')
	filelist.sort()
	dat = []
	for filename in filelist: dat.append(angles_calcul(filename[:-4]))
	return dat

#Defining the arguments:
parser = argparse.ArgumentParser(description="Get the secondary structure from the phi and psi angles.")
parser.add_argument("files", help="The .npy file with the angles or a folder with the PDB files to calculate the angles.")
parser.add_argument("-alpha", "-a", action = "store_true", default = False,
					help="When present  the alpha  helix region will be studied.")
parser.add_argument("-beta", "-b", action = "store_true", default = False,
					help="When present  the beta strand  region will be studied.")
parser.add_argument("-polyproline", "-ppii", action = "store_true", default = False,
					help="When present  the polyproline II helix  region will be studied..")
parser.add_argument("-stride", "-st", nargs = '+', default = False,
					help = "The  PDB files directory to calculate the stride, followed by the desired structure: 'alpha' and/or 'beta'.")
parser.add_argument("-save_figure","-sf", nargs = "+", default = False,
					help = "The path where the figure/s will be saved (including  prefix, see documentation for more details).")
parser.add_argument("-save_numpy","-sn", nargs = "+", default = False,
					help = "The path to save the numpy/s array.")
parser.add_argument("-structure_definition", "-sd", choices=["profasi","pappu","blackledge"],
					default="blackledge",
					help = "Which regions definition you want to use: ''profasi, 'pappu' or 'blackledge'. by default 'blackledge'.\nYou can define you own region with the customizer_region option.")
parser.add_argument("-groups","-g", type = int, nargs = "+", default = False,
					help = "The initial and the final groups to show. The default values are all groups except the group zero (which indicates the aminoacids without the desired conformation).")
parser.add_argument("-residues","-r", type = int, nargs = "+", default = False,
					help = "The initial and the final residues to show. By default all the residues are shown excep the first one and the last one (for more details see documentation).")
parser.add_argument("-hr","-helix-per-residue", action = "store_true", default = False,
					help = "When present the program draws the alpha-helix percentage per residue.")
parser.add_argument("-hgt","-helix-per-group", action = "store_true", default = False,
					help = "When present the program draws the alpha helix region length per residue and temperature.")
parser.add_argument("-txt", default = False,
					help = "When present the program will save the percentages in a .txt file. It should indicate the path to save the file.")
parser.add_argument("-temp","-temperature", default = False, nargs = '+',
					help = "This option set the temperatures shown in the y axe.")
parser.add_argument("-rgc","-range_colorbar", default = False, nargs = 2,
					help = "This option sets the minimum and maximum percentage shown in the resulting figure.")
parser.add_argument("-w","-weights", default = False, nargs = '+',
					help = "This option specifies the complete path to the weights file.")
parser.add_argument("-cm","-color_map", default = "jet", choices = ['jet','binary'],
					help = "This option specifies the colormap to plot the images for black and white the option should be: binary.")
parser.add_argument("-customized_region", "-cr", default = False, nargs = 5,
					help = "This option defines a customized region in the Ramachandran Plot.\nUsage = -cr  Structure  phi0 phi1 psi0 psi1.\nWhere Structure is the conformation's name in the region  and phi/psi0 is the minimum value and the phi/psi1 is the maximun value for the angle.")

global args
args = parser.parse_args()

if not figures: args.txt = True

if args.polyproline and args.structure_definition == "profasi":
	args.structure_definition = "blackledg"

global all_data
if args.files.split(".")[-1] == "npy":
	all_data = np.load(args.files)
elif args.files.split("/")[-1] == "":
	all_data = np.asarray(pdb_npy(args.files))
else:
	print "Incorrect data type.\nThis program only takes as valid input a numpy array or a folder with multple PDB files."
	sys.exit()

global weights
if args.w:
	type = args.w.split(".")[-1]
	if type == "txt":
		weights = np.loadtxt(args.w)
	elif type == "npy":
		weights = np.load(args.w)
	else:
		print "This program only takes the weights from a .txt file or a .npy file"
	if temporary.shape[0] != all_data.shape[0]:
		print "The weights and the number of structures do not match."
		sys.exit()
	weights /= weights.sum()
else:
	weights = np.ones(all_data.shape[0])

if not args.residues: args.residues = [0,all_data.shape[1]]
else: args.residues = [args.residues[0]-2,args.residues[1]-1]
if not args.groups: args.groups = [1,all_data.shape[1]+1]

aminoacids = all_data.shape[1]
global all_structure
if args.structure_definition == "blackledge"  : 
	all_structure = np.asarray([blackledge(data) for data in all_data])
elif args.structure_definition == "profasi":
	all_structure = np.asarray([profasi(data) for data in all_data])
elif args.structure_definition == "pappu":
	all_structure = np.asarray([pappu(data) for data in all_data])
elif args.customized_region:
	all_structure = np.asarray([custom(data) for data in all_data])

if args.stride:
	"""Calculating stride"""
	pdb_list = glob.glob(args.stride[0]+"*.pdb")
	if not pdb_list: print "Directory no valid: "+args.stride[0]
	else:
		print "Calculating stride with the pdbs in directory: %s"%args.stride[0]
		all_structure_stride = np.asarray([stride(pdb[:-4]) for pdb in pdb_list])
		print "Stride completed"
		if args.stride[1] == 'alpha' or args.stride[2] == 'alpha':
			d_stride_alpha = np.zeros([all_structure_stride.shape[1], all_structure_stride.shape[1]+1])
			for structure in all_structure_stride: d_stride_alpha = d_stride_alpha + count(structure,"alpha")
			stride_alpha_percentage = d_stride_alpha/all_structure_stride.shape[0]
			if figures:
				images(stride_alpha_percentage, "-stride-alpha")
			if args.save_numpy:
				np.save(args.save_numpy[0]+"-stride-beta-strand-percentage", d_stride[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]])
		elif args.stride[1] == 'beta' or args.stride[2] == 'beta':
			d_stride_beta = np.zeros([all_structure_stride.shape[1], all_structure_stride.shape[1]+1])
			for structure in all_structure_stride: d_stride_beta = d_stride_beta + count(structure,"beta")
			stride_beta_percentage = d_stride_beta/all_structure_stride.shape[0]
			if figures:
				images(stride_alpha_percentage, "-stride-beta")
			if args.save_numpy:
				np.save(args.save_numpy[0]+"-stride-alpha-helix-percentage", d_stride[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]])

if args.alpha:
	d_alpha = np.zeros([aminoacids, aminoacids+1])
	for structure,w in zip(all_structure,weights):
		d_alpha += w*count(structure,"alpha")
	alpha_percentage = d_alpha / all_structure.shape[0]
	if figures:
		images(alpha_percentage, "alpha")
	if args.save_numpy:
		np.save(args.save_numpy[0]+"-alpha-helix-percentage-%s-definition"%args.structure_definition, d_alpha[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]])

if args.beta:
	d_beta = np.zeros([aminoacids, aminoacids+1])
	for structure,w in zip(all_structure,weights):
		d_beta += w*count(structure,"beta")
	beta_percentage = d_beta / all_structure.shape[0]
	if figures:
		images(beta_percentage, "beta")
	if args.save_numpy:
		np.save(args.save_numpy[0]+"-beta-strand-percentage-%s-definition"%args.structure_definition, d_beta[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]])

if args.polyproline:
	d_ppii = np.zeros([aminoacids, aminoacids+1])
	for structure,w in zip(all_structure,weights):
		d_ppii += w * count(structure,"PPII")
	ppii_percentage = d_ppii / all_structure.shape[0]
	if figures:
		images(ppii_percentage, "polyprolineII")
	if args.save_numpy:
		np.save(args.save_numpy[0]+"polyproline-percentage", d_ppii[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]])

if args.customized_region:
	d_custom = np.zeros([aminoacids, aminoacids+1])
	for structure,w in zip(all_structure,weights):
		d_custom += w * count(structure,args.customized_region)
	custom_percentage = d_custom / all_structure.shape[0]
	if figures:
		images(custom_percentage, "custom region")
	if args.save_numpy:
		np.save(args.save_numpy[0]+"-custom-conformation-percentage-%s-definition"%args.structure_definition, d_custom[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]])

if args.hr and figures:
	try:
		d_alpha
	except NameError:
		d_alpha = np.zeros([aminoacids, aminoacids+1])
		for structure,w in zip(all_structure,weights): d_alpha = d_alpha + w*count(structure,"alpha")
		d_alpha = d_alpha/all_structure.shape[0]
	pl.figure()
	pl.plot(d_alpha[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]].sum(axis=1))
	#pl.title("% Alpha-helix", fontsize =12)
	pl.xlim(1,d_alpha.shape[0])
	if args.save_figure:
		pl.savefig(args.save[0]+"-alpha-helix-percentage-per-residue-%s-definition"%args.structure_definition)
	if args.save_numpy:
		np.save(args.save_numpy[0]+"-alpha-helix-percentage-per-residue-%s-definition"%args.structure_definition, d_alpha[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]].sum(axis=1))

if args.hgt:
	all_files = glob.glob("-".join(args.files.split('-')[:-1])+"*.npy")
	#print "All the files should be from the same type: .npy."
	all_files.sort()
	all_data_hgt = []
	all_data_hgt = np.asarray([np.load(file) for file in all_files])
	global all_structure_hgt
	all_structure_hgt = []
	for data in all_data_hgt:
		if args.structure_definition == "blackledge"  : 
			all_structure_hgt.append(np.asarray([blackledge(dat) for dat in data]))
		elif args.structure_definition == "profasi":
			all_structure_hgt.append(np.asarray([profasi(dat) for dat in data]))
		elif args.structure_definition == "pappu":
			all_structure_hgt.append(np.asarray([pappu(dat) for dat in data]))
	all_structure_hgt = np.asarray(all_structure_hgt)
	d0 = np.zeros([aminoacids, aminoacids+1])
	for struct in all_structure_hgt[0]: d0 += count(struct,"alpha")
	d0 /= all_structure_hgt[0].shape[0]
	d = d0.sum(axis=0)/d0.shape[1]
	for structure in all_structure_hgt[1:]:
		d0 = np.zeros([aminoacids, aminoacids+1])
		for struct in structure: d0 += count(struct,"alpha")
		d0 /= structure.shape[0]
		d = np.row_stack((d,d0.sum(axis=0)/d0.shape[1]))
	if figures:
		temperatures = []
		if not args.temp:
			for i in range(len(all_files)): temperatures.append(i)
		else:
			for element in args.temp: temperatures.append(element)
		fig = pl.figure("Structured region lenght per residue and temperature")
		ax = fig.add_subplot(111)
		cs = ax.contourf(range(args.groups[0],args.groups[1]),np.asarray(temperatures),
						d[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]],
						np.linspace(0, d[:,1:].max(), 1000))
		if args.rgc: cs.set_clim(args.rgc[0],args.rgc[1])
		cb=fig.colorbar(cs)
		cb.set_ticks(np.linspace(0,np.round(d[:,1:].max(),decimals = 2), 11))
		(ydim, xdim) =d[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]].shape
		starting_residue=2
		xticks_spacing=int(ax.get_xticks()[1]-ax.get_xticks()[0])
		yticks_spacing=int(ax.get_xticks()[1]-ax.get_yticks()[0])
		ax.grid()
		fig.show()
	if args.save_figure:
		pl.savefig("%s-alpha-helix-percentage-per-group-and-temperature-%s"%(args.save_figure[0],args.structure_definition))
	if args.save_numpy:
		np.save(args.save_numpy[0]+'structured-region-lenght-with-temperature', d[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]])

if figures: pl.show()

if args.txt and (args.stride or args.alpha or args.beta or args.polyproline or args.hr or args.hrt or args.hgt):
	if args.stride:
		try:
			stride_alpha_percentage
			np.savetxt(args.txt+"ss-map-stride-alpha-percentage-%s-definition.txt"%args.structure_definition,stride_alpha_percentage,fmt = '%4f')
		except NameError: pass
		try:
			stride_beta_percentage
			np.savetxt(args.txt+"ss-map-stride-beta-percentage-%s-definition.txt"%args.structure_definition,stride_beta_percentage,fmt = '%4f')
		except NameError: pass
	if args.alpha:
		np.savetxt(args.txt+"ss-map-alpha-percentage%s-definition.txt"%args.structure_definition,alpha_percentage,fmt = '%4f')
	if args.beta:
		np.savetxt(args.txt+"ss-map-beta-percentage%s-definition.txt"%args.structure_definition,beta_percentage,fmt = '%4f')
	if args.polyproline:
		np.savetxt(args.txt+"ss-map-polyproline-percentage%s-definition.txt"%args.structure_definition,ppii_percentage,fmt = '%4f')
	if args.hr:
		d_hr = alpha_percentage[args.residues[0]:args.residues[1], args.groups[0]:args.groups[1]].sum(axis=1)
		np.savetxt(args.txt+"ss-map-helix-per-residue%s-definition.txt"%args.structure_definition,d_hr,fmt = '%4f')
elif args.txt and not (args.stride or args.alpha or args.beta or args.polyproline or args.hr or args.hrt or args.hgt):
	 print "You have calculated nothing."

