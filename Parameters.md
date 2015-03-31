# Parameters #



## Mandatory parameter ##

This program has only one mandatory parameter which is the path to either the numpy array with the  φ and ψ angles, or to the folder containing the PDB files.

## Optional parameters ##

### Accepted conformations ###

The following options set the conformations to look for in the ensemble.
  * **-alpha** or **-a**: When present this option makes the program study the alpha helix region.
  * **-beta** or **-b**: When present this option makes the program study the beta strand helix region.
  * **-polyproline** or **-ppii**:  When present this option makes the program study the polyproline II helix region.

### Ramachandran plot regions ###

The regions in the Ramachandran plot are defined by the following options:
  * **-structure\_definition** or **-sd** _{profasi, pappu, blackledge}_: This option sets the divisions using on of the following predefined arguments:
    * _blackledge_: this option sets the regions defined in the article **_Mapping the Potential Energy Landscape of Intrinsically Disordered Proteins at Amino Acid Resolution_**, Valery Ozenne, Robert Schneider, Mingxi Yao, Jie-rong Huang, Loic Salmon, Markus Zweckstetter,Malene Ringkjobing Jensen, and Martin Blackledge, _**J.Am.Chem.Soc**, 2012, 134, 15138-15148_ in Figure 2 page 15140.
    * _pappu_: this option sets the regions defined in the article **_Net charge per residue modulates conformational ensembles of intrinsically disordered proteins_**, Albert H. Mao, Scott L. Crick, Andreas Vitalis, Caitlin L. Chicoine, and Rohit V. Pappu, _**PNAS**, 2010, vol. 107, no 18, 8183-8188_ at the section: Supporting information image Secondary structure definitons.
  * _profasi_: this option sets the default regions used by the program PROFASI
    * **-customized\_region** or **-cr** _CONFORMATION PHI0 PHI1 PSI0 PSI1_: This option allows the user to define a rectangular region in the Ramachandran plot and assign it any conformation the user wants. It takes five (5) arguments: conformation phi0 phi1 psi0 psi1. The argument _conformation_ is the region's conformation name _phi/psi0_ is the minimum value and _phi/psi1_ is the maximum value for the angles.

### STRIDE program ###

The ss-map program can call the program STRIDE to calculate the conformations with the option:
  * **-stride** or **-st** _PDB\_DIRECTORY {alpha, beta}_: When using this option the only accepted conformations are _alpha_ or _beta_. This options takes a minimum of two arguments and a maximum of three. The first argument must be the directory with the PDB files to call the STRIDE, the second and third arguments should be _alpha_, _beta_ or both dependign on which structure you want to study.

### Images properties ###

The ss-map allows the user to set some properties of the image with the following options:
  * **-color\_map** or **-cm** _{jet, binary}_: This option allows the user to choose between two (2) arguments _jet_ and _binary_, the first one is the default, the second one allows the user to create the images in black and white.
  * **-range\_colorbar** or **-rgc** _MINIMUM MAXIMUM_: This option sets the colorbar range. It takes two (2) arguments the _minimum_ and the _maximum_ percentage to be used in the resulting figure/s.

### Data properties ###

The following options allow the user to set the data to show/save, it affects the figures and the data to save.
  * **-groups** or **-g** _MINIMUM\_LENGHT MAXIMUM LENGHT_: This option sets the minimum and the maximum lenght to show/save. By default this option is set to show from the single aminoacids (one) to the protein's lenght. The minimum length accepted is zero which represents the aminoacids without the desired conformation.
  * **-residues** or **-r** _FIRST\_RESIDUE LAST\_RESIDUE_: This option sets the first and the last residue to show/save. By default this options is set to show all de aminoacids in the protein, except the first and the last, for more details on why these aminoacids cannot be shown see the [FAQs](FAQs.md)


### Weightening the ensemble ###

This option allows the user to assign a weigth to every structure in the ensemble.
  * **-weigths** or **-w** _WEIGTHS\_FILE_: This option indicates the path to the file with the ensemble weights. The accepted files are either a numpy array with extension .npy or .npz, or a text file with extension .txt.

### Other Plots ###

The following options allow the user to create another type of images. The following options accept the different [Ramachandran plot regions](Parameters#Ramachandran_plot_regions.md), [Images](Parameters#Images_properties.md) and [Data](Parameters#IData_properties.md) properties, but only work with the alpha conformation.
  * **-helix-per-residue** or **-hr**: This option generates an image which represents the percentage of alpha-helix per residue.
  * **-helix-per-group** or **-hgt**: This option allows the user to represent in a single image different ensembles. To use this option the angles φ and ψ must be already calculated and saved in a numpy array for each ensemble. All the arrays must have the same name except for a sufix separated by the simbol '-', i.e.: arraysensemble-1.npy, ..., arraysensemble-n.npy .
  * **-temperature** or **-temp** [TEMPERATURE\_2 ... TEMPERATURE\_n](TEMPERATURE_1.md): This option specifies the temperatures to show on the axe y in the figure generated by the **-hgt** option.


### Saving the results ###

The following options allow the user to automatically save the figures and/or the data.
  * **-save\_figure** or **-sf** _PATH\_TO\_SAVE_: This option specifies the path and the prefix to save the figure/s. For more details on how to use this option see [Saving the figures](Usage#Saving_the_figures.md)
  * **-save\_numpy** or **-sn** _PATH\_TO\_SAVE_: This option specifies the path and the prefix to save the numpy array/s with the percentage information. One array per conformation will be saved, for more details see [Saving the arrays](Usage#Saving_the_arrays.md)
  * **-txt** _PATH\_TO\_CREATE\_THE\_TXT\_FILE/S_: This option specifies the path to save the percentages in a .txt file. The files name will be set automatically using the following structure: ss-map-stride(if used)-conformation(_alpha, beta_ or_polyprolyne_)-structure\_definition (_profasi, pappu, blackledge_ or _customized_)-definition.txt

### Help ###
The program also has a minor help implemented.
  * **--help** or **-h**: Shows the program help message.