#!/usr/bin/env python3
#
#============================ SS-map ==================================
# ss-map is a Python program to visualize the proteins ensembles
# secondary structure.
# ss-map is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# ss-map is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with ss-map.  If not, see <http://www.gnu.org/licenses/>.
# Written by: Jelisa Iglesias.
#
#======================================================================
#
# Authors:
# Jelisa Iglesias
# Ramon Crehuet
#======================================================================
# For more information visit:
# https://github.com/rcrehuet/ss-map
#
# This work is done by the Theoretical and Computational Group
# of the IQAC (CSIC)
# http://iqac.csic.es/qteor
#


"""

The numpy library is required. But he pylab library is
optional, if it is not present no images will be generated
and the information will be saved in a .txt file
"""
figures = True
import sys
import os
import argparse
import subprocess as subp
import glob
import numpy as np

try: 
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter #Needed to zoom in figure
except ImportError:
    print("You do not haveinstalled the Matplotlib module.\n\
    The program will not generate any image, if you want to save the data \
    use any of the following options:\n-save_numpy \n-txt")
    figures = False
"""
The following libraries are installed when intalling python as part
of the Standard Library
"""

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
        if element[0] >= 0:
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

def count(data, struct):
    """
    This function counts how many aminoacids are in a given structure in a row.
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
    def xformat(x, pos):
        """
        Format the length of the fragment, as an integer starting by 1
        """
        return '{:.0f}'.format(x+starting_group,)

    def yformat(y, pos):
        """
        Format the the residue number
        """    
        return '{:.0f}'.format(y+starting_residue,)
    
    xformatter = FuncFormatter(xformat)
    yformatter = FuncFormatter(yformat)
    fig = plt.figure(structure)
    ax = fig.add_subplot(111)
    cs = ax.matshow(percentages[args.residues[0]:args.residues[1], \
         args.length[0]:args.length[1]], cmap = args.cm)
    if args.rgc: cs.set_clim(float(args.rgc[0]),float(args.rgc[1]))
    fig.colorbar(cs)
    (ydim, xdim) =percentages[args.residues[0]:args.residues[1], \
                  args.length[0]:args.length[1]].shape
    starting_residue=args.residues[0]+2
    starting_group = args.length[0]
    ax.yaxis.set_major_formatter(yformatter)
    ax.xaxis.set_major_formatter(xformatter)
    ax.grid()
    fig.show()
    if args.save_figure:
        plt.savefig("%s-%s-estructure-%s-definition"%(args.save_figure,\
        structure,args.structure_definition))
    return

def numpys (structure, results):
    np.save("%s-%s-percentage-%s-definition"%(args.save_numpy, \
    structure,args.structure_definition), \
    results[args.residues[0]:args.residues[1], args.length[0]:args.length[1]])
        
#Defining the arguments:
parser = argparse.ArgumentParser(description="Get the secondary structure from the phi and psi angles.")
parser.add_argument("files", 
    help="The .npy file with the angles (calculated with get_angles.")
conformations = parser.add_argument_group("conformations", "All the predefined accepted conformations.")
conformations.add_argument("--alpha", "-a", action = "store_true", default = False,
    help="When present  the alpha  helix region will be studied.")
conformations.add_argument("--beta", "-b", action = "store_true", default = False,
    help="When present  the beta strand  region will be studied.")
conformations.add_argument("--polyproline", "-ppii", action = "store_true", default = False,
    help="When present  the polyproline II helix  region will be studied..")
ramachandran_regions = parser.add_argument_group("Ramachandran regions", "The commands to set the ramachandran regions.")
ramachandran_regions.add_argument("--structure_definition", "-sd", choices=["profasi","pappu","blackledge"],
                    default="blackledge",
    help = "Which regions definition you want to use: ''profasi, 'pappu' or 'blackledge'. by default 'blackledge'.\nYou can define you own region with the customizer_region option.")
ramachandran_regions.add_argument("--customized_region", "-cr", default = False, nargs = 5,
    help = "This option defines a customized region in the Ramachandran Plot. Usage = -cr  Structure  phi0 phi1 psi0 psi1. Where Structure is the conformation's name in the region  and phi/psi0 is the minimum value and the phi/psi1 is the maximun value for the angles.")

stride2 = parser.add_argument_group("STRIDE", "The command that calls the external program STRIDE.")
stride2.add_argument("--stride", "-st", default = False,
    help = "The path to STRIDE executable.")

images_properties = parser.add_argument_group("Images properties","The commands to change the images properties.")
images_properties.add_argument("--color_map", "-cm", dest='cm', default = "Blues",
    help = "This option specifies the colormap to plot the images. \
    For black and white the option should be: binary.")
images_properties.add_argument("--range_colorbar", "-rgc", dest='rgc', default = False, nargs = 2,
    help = "This option sets the minimum and maximum percentage shown in the resulting image.")
data_properties = parser.add_argument_group("Data properties","The commands to change the data shown.")
data_properties.add_argument("--length","-l", type = int, nargs = 2, default = False,
    help = "The initial and the final lengths to show. The default values are all lengths except the zero (which indicates the aminoacids without the desired conformation).")
data_properties.add_argument("--residues","-r", type = int, nargs = 2, default = False,
    help = "The initial and the final residues to show. By default all the residues are shown.")

ensemble_weights = parser.add_argument_group("Ensemble weights","The commands to specify the ensembles weights.")
ensemble_weights.add_argument("--weights","-w", dest='w',default = False,
    help = "This option specifies the complete path to the weights file.")

other_plots = parser.add_argument_group("Other plots", "The commands to generate other plots.")
other_plots.add_argument("--helix-per-residue","-hr", dest='hr', action = "store_true", default = False,
    help = "When present the program draws the alpha-helix percentage per residue.")
other_plots.add_argument("--helix-per-length", "-hl", dest='hg', action = "store_true", default = False,
    help = "When present the program draws the alpha helix region length.")
#other_plots.add_argument("-temp","-temperature", default = False, nargs = '+',
#    help = "This option set the temperatures shown in the y axe.")

saving_the_results = parser.add_argument_group("Saving the results","Commands to save the data.")
saving_the_results.add_argument("--save_figure","-sf",  default = False,
    help = "The path where the figure/s will be saved (including  prefix, see documentation for more details).")
saving_the_results.add_argument("--save_numpy","-sn",  default = False,
    help = "The path to save the numpy/s array.")
saving_the_results.add_argument("-txt", default = False,
    help = "When present the program will save the percentages in a .txt file. It should indicate the path to save the files.")
saving_the_results.add_argument("-no_fig",action = "store_true", default = False,
    help = "When present the program will not generate figures (only meaningful when you save data.")


args = parser.parse_args()

if args.no_fig : figures = False
#if not figures: args.txt = True

if args.polyproline and args.structure_definition == "profasi":
    args.structure_definition = "blackledg"

if args.files.split(".")[-1] == "npy":
    all_data = np.load(args.files)
else:
    print("Incorrect data type.\nThis program only takes as valid input a numpy array.")
    sys.exit()

if args.w:
    typo = args.w.split(".")[-1]
    if typo == "txt":
        weights = np.loadtxt(args.w)
    elif typo == "npy":
        weights = np.load(args.w)
    else:
        print("This program only takes the weights from a .txt file or a .npy file.")
    if weights.shape[0] != all_data.shape[0]:
        print("The number of weights and the number of structures do not match.")
        sys.exit()
    weights /= weights.sum()
    weights *= len(weights)
else:
    weights = np.ones(all_data.shape[0])

if not args.residues: args.residues = [0,all_data.shape[1]]
else: args.residues = [args.residues[0]-2,args.residues[1]-2]
if args.length:
    args.length[-1] += 1
else: 
    args.length = [1,all_data.shape[1]+1]

aminoacids = all_data.shape[1]
if args.structure_definition == "blackledge"  : 
    all_structure = np.asarray([blackledge(data) for data in all_data])
elif args.structure_definition == "profasi":
    all_structure = np.asarray([profasi(data) for data in all_data])
elif args.structure_definition == "pappu":
    all_structure = np.asarray([pappu(data) for data in all_data])
if args.customized_region:
    all_structure = np.asarray([custom(data) for data in all_data])
    args.structure_definition = "customized"

if args.stride:
    """Calculating stride"""
    pdb_list = glob.glob(args.stride[0]+"*.pdb")
    if not pdb_list: print("Directory no valid: "+args.stride[0])
    else:
        alfa = False
        beta = False
        print("Calculating stride with the pdbs in directory: %s"%args.stride[0])
        all_structure_stride = np.asarray([stride(pdb[:-4]) for pdb in pdb_list])
        print("Stride completed")
        try: args.stride[1] == 'alpha' or args.stride[2] == 'alpha'
        except IndexError: pass
        else: alfa = True
        if alfa:
            d_stride_alpha = np.zeros([all_structure_stride.shape[1], all_structure_stride.shape[1]+1])
            for structure in all_structure_stride: d_stride_alpha = d_stride_alpha + count(structure,"alpha")
            stride_alpha_percentage = d_stride_alpha/all_structure_stride.shape[0]
            if figures:
                images(stride_alpha_percentage, "stride-alpha-helix")
            if args.save_numpy:
                numpys("stride-alpha-helix", stride_alpha_percentage)
        try: args.stride[1] == 'beta' or args.stride[2] == 'beta'
        except IndexError: pass
        else: beta = True
        if beta:
            d_stride_beta = np.zeros([all_structure_stride.shape[1], all_structure_stride.shape[1]+1])
            for structure in all_structure_stride: d_stride_beta = d_stride_beta + count(structure,"beta")
            stride_beta_percentage = d_stride_beta/all_structure_stride.shape[0]
            if figures:
                images(stride_beta_percentage, "stride-beta-strand")
            if args.save_numpy:
                numpys("stride-beta-strand", stride_beta_percentage)

if args.alpha:
    d_alpha = np.zeros([aminoacids, aminoacids+1])
    for structure,w in zip(all_structure,weights):
        d_alpha += w*count(structure,"alpha")
    alpha_percentage = d_alpha / all_structure.shape[0]
    if figures:
        images(alpha_percentage, "alpha-helix")
    if args.save_numpy:
        numpys("alpha-helix", alpha_percentage)

if args.beta:
    d_beta = np.zeros([aminoacids, aminoacids+1])
    for structure,w in zip(all_structure,weights):
        d_beta += w*count(structure,"beta")
    beta_percentage = d_beta / all_structure.shape[0]
    if figures:
        images(beta_percentage, "beta-strand")
    if args.save_numpy:
        numpys("beta-strand", beta_percentage)

if args.polyproline:
    d_ppii = np.zeros([aminoacids, aminoacids+1])
    for structure,w in zip(all_structure,weights):
        d_ppii += w * count(structure,"PPII")
    ppii_percentage = d_ppii / all_structure.shape[0]
    if figures:
        images(ppii_percentage, "polyprolineII-helix")
    if args.save_numpy:
        numpys("polyprolineII-helix", ppii_percentage)

if args.customized_region:
    d_custom = np.zeros([aminoacids, aminoacids+1])
    for structure,w in zip(all_structure,weights):
        d_custom += w * count(structure,args.customized_region)
    custom_percentage = d_custom / all_structure.shape[0]
    if figures:
        images(custom_percentage, "custom-region")
    if args.save_numpy:
        numpys("custom-region", custom_percentage)

if args.hr:
    try:
        d_alpha
    except NameError:
        d_alpha = np.zeros([aminoacids, aminoacids+1])
        for structure,w in zip(all_structure,weights): d_alpha = d_alpha + w*count(structure,"alpha")
        alpha_percentage = d_alpha/all_structure.shape[0]
    if figures:
        plt.figure()
        plt.plot(alpha_percentage[args.residues[0]:args.residues[1], args.length[0]:args.length[1]].sum(axis=1))
    #plt.title("% Alpha-helix", fontsize =12)
        plt.xlim(1,alpha_percentage.shape[0])
    if args.save_figure:
        plt.savefig(args.save_figure+"-per-residue-%s"%args.structure_definition)
    if args.save_numpy:
        np.save(args.save_numpy+\
        "-per-residue-%s"\
        %args.structure_definition,\
        alpha_percentage[args.residues[0]:args.residues[1], \
        args.length[0]:args.length[1]].sum(axis=1))

if args.hg:
    try:
        d_alpha
    except NameError:
        d_alpha = np.zeros([aminoacids, aminoacids+1])
        for structure,w in zip(all_structure,weights): d_alpha = d_alpha + w*count(structure,"alpha")
        alpha_percentage = d_alpha/all_structure.shape[0]
    if figures:
        plt.figure()
        plt.plot(alpha_percentage[args.residues[0]:args.residues[1], args.length[0]:args.length[1]].sum(axis=0))
        #plt.title("% Alpha-helix", fontsize =12)
        plt.xlim(1,alpha_percentage.shape[1])
    if args.save_figure:
        plt.savefig(args.save_figure+"-per-length-%s"%args.structure_definition)
    if args.save_numpy:
        np.save(args.save_numpy+\
        "-per-length-%s"\
        %args.structure_definition,\
        alpha_percentage[args.residues[0]:args.residues[1], \
        args.length[0]:args.length[1]].sum(axis=0))
        
if figures: plt.show()
#else: args.txt = True

if args.txt and (args.stride or args.alpha or args.beta or args.polyproline or args.hr or args.hrt or args.hgt):
    if args.stride:
        try:
            stride_alpha_percentage
            np.savetxt(args.txt+"ss-map-stride-alpha-percentage-%s-definition.txt"\
            %args.structure_definition,stride_alpha_percentage,fmt = '%4f')
        except NameError: pass
        try:
            stride_beta_percentage
            np.savetxt(args.txt+"ss-map-stride-beta-percentage-%s-definition.txt"\
            %args.structure_definition,stride_beta_percentage,fmt = '%4f')
        except NameError: pass
    if args.alpha:
        np.savetxt(args.txt+"ss-map-alpha-percentage%s-definition.txt"\
        %args.structure_definition,alpha_percentage,fmt = '%4f')
    if args.beta:
        np.savetxt(args.txt+"ss-map-beta-percentage%s-definition.txt"\
        %args.structure_definition,beta_percentage,fmt = '%4f')
    if args.polyproline:
        np.savetxt(args.txt+"ss-map-polyproline-percentage%s-definition.txt"\
        %args.structure_definition,ppii_percentage,fmt = '%4f')
    if args.hr:
        d_hr = alpha_percentage[args.residues[0]:args.residues[1], args.length[0]:args.length[1]].sum(axis=1)
        np.savetxt(args.txt+"ss-map-helix-per-residue%s-definition.txt"\
        %args.structure_definition,d_hr,fmt = '%4f')
elif args.txt and not (args.stride or args.alpha or args.beta or args.polyproline \
   or args.hr or args.hrt or args.hgt):
     print("You have calculated nothing.")

