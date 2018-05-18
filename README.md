# SS-map

This project aims to develope a program to visualize the proteins ensembles secondary structure in a simple way.

The ss-map is a Python program to visualize the proteins ensembles secondary structure, using the angles φ and ψ. The program takes either a folder with multiple PDB files (one for each protein in the ensemble) with which calculates the angles φ and ψ, or an array with the angles φ and ψ for each structure in the ensemble. It returns either an image, a numpy array and/or a .txt file containing a matrix (or graphical a representation of this matrix) which shows in how many structures of the ensemble (in %) the residue y is forming a structured region of lenght x.

## Requirements:

    Python language with libraries:
        matplotlib (only necessary to get the images)
        mdtraj (if trajectories or pdb files need to be read)
        Numpy 

Important: If the library Matplotlib isn't installed the program will write the information in a .txt file but will not produce any images.

## Citation

The science behind these code and its use is described in a publication submitted to Intrinsically Disordered Proteins journal. If you use SS-map, please cite this publication (when of if accepted, details will be updated):

[Jelisa Iglesias, Melchor Sanchez-Martinez, Ramon Crehuet 'SS-map: Visualizing cooperative secondary structure elements in protein ensembles', *Intrinsically Disordered Proteins* 2013; 1.](http://www.tandfonline.com/doi/full/10.4161/idp.25323) 
