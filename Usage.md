# Usage #



This section explains how to use the SS-map program different options. The SS-map program can be used either in Python or in Windows, whenever the python language is installed. The program can be executed using one of the following instructions in a console.
```
$ ./ss-map.py the/path/to/the/numpy_array.npy     (1)
or
$ ./ss-map.py the/path/to/the/pdb/files/          (2)
```
To simplify the text, throughout this wiki the first instruction (1) will be used.

Just with one of this instructions the program will do nothing. To generate any data the user should specify at least one of the accepted conformations to study, or the option _-helix\_per\_group_.


## Choosing the conformation ##

The SS-map can study three (3) different conformations: _alpha-helix, beta-strand_and_polyprolineII helix_.
To generate an image representing the alpha-helix conformation the instruction should be:
```
$ ./ss-map.py the/path/to/the/numpy_array.npy -alpha         (4)
or 
$ ./ss-map.py the/path/to/the/numpy_array.npy -a             (5)
```
Both instructions generate the images in **Figure 1**
| ![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-alpha-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-alpha-blackledge-definition.png) | ![https://ss-map.googlecode.com/svn/wiki/images/examples/beta-hairpin-estructura-alpha-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/beta-hairpin-estructura-alpha-blackledge-definition.png) |
|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 1**: Both images have been generated with the instruction above and represent the alpha-helix content of two different ensembles, the left one represents the ensemble alpha\_helix\_angles\_array.npy and the right one the ensemble beta\_hairpin\_angles\_array.npy (both downloadables from [here](https://drive.google.com/folderview?id=0B5Wxvh8ya41STDlra2pGZkUxaEk&usp=sharing))


To generate an image representing the beta-strand conformation the instruction should be:
```
$ ./ss-map.py the/path/to/the/numpy_array.npy -beta          (6)
or 
$ ./ss-map.py the/path/to/the/numpy_array.npy -b             (7)
```
Both instructions generate the images in **Figure 2**
| ![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-beta-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-beta-blackledge-definition.png) | ![https://ss-map.googlecode.com/svn/wiki/images/examples/beta-hairpin-estructura-beta-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/beta-hairpin-estructura-beta-blackledge-definition.png) |
|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 2**: Both images have been generated with the instruction above and represent the beta-strand content of two different ensembles, the left one represents the ensemble alpha\_helix\_angles\_array.npy  and the right one the ensemble beta\_hairpin\_angles\_array.npy (both downloadables from [here](https://drive.google.com/folderview?id=0B5Wxvh8ya41STDlra2pGZkUxaEk&usp=sharing))

To generate an image representing the polyprolineII helix conformation the instruction should be:
```
$ ./ss-map.py the/path/to/the/numpy_array.npy -polyproline   (8)
or 
$ ./ss-map.py the/path/to/the/numpy_array.npy -ppii          (9)
```
Both instructions generate the images in **Figure 3**
| ![https://ss-map.googlecode.com/svn/wiki/images/examples/Polyproline.png](https://ss-map.googlecode.com/svn/wiki/images/examples/Polyproline.png) |
|:--------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 3**

All this options can be used at the same time over one ensemble and the program will produce the same number of images as the number of conformations specified, at the same time.
```
$ ./ss-map.py the/path/to/the/numpy_array.npy -a - beta      (10)
```

| ![https://ss-map.googlecode.com/svn/wiki/images/examples/beta-hairpin-both-estructura-alpha-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/beta-hairpin-both-estructura-alpha-blackledge-definition.png)| ![https://ss-map.googlecode.com/svn/wiki/images/examples/beta-hairpin-both-estructura-beta-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/beta-hairpin-both-estructura-beta-blackledge-definition.png)|
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 4**: Both images were generated at the same time, with just one instruction.


## Using the different regions ##

The default regions are the ones defined by the option Blackledge (for more information on the regions see [Ramachandran plot regions](Parameters#Ramachandran_plot_regions.md)). To change this regions the user has two options:
  1. To use one of the predefined divisions by the program using the option **structure\_definition**:
```
$ ./ss-map.py the/path/to/the/numpy_array.npy -a -structure_definition blackledge        (11)
or 
$ python ss-map the/path/to/the/numpy_array.npy -a -sd blackledge                        (12)

$ ./ss-map.py the/path/to/the/numpy_array.npy -a                                         (13)

# The -sd option accepts three options: blackledge, pappu and profasi. And each 
# of this options defines a different partition of the Ramachandran Space (for more details see the Ramachandran plot
# regions section in Parameters).

$ ./ss-map.py the/path/to/the/numpy_array.npy -a -sd pappu                               (14)

$ ./ss-map.py the/path/to/the/numpy_array.npy -a -sd profasi                             (15)

```
  1. To define a rectangular area and assign it a conformation using the option **customized\_region**:
```
$ ./ss-map.py the/path/to/the/numpy_array.npy -customized_region conf -175 -30 -70 16    (16)
or 
$ ./ss-map.py the/path/to/the/numpy_array.npy -cr conf -175 -30 -70 160                  (17)
```

| ![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-alpha-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-alpha-blackledge-definition.png) |
|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 5**: All three instructions (11), (12) and (13) generate this image.

| ![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-alpha-pappu-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-alpha-pappu-definition.png) | ![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-alpha-profasi-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-alpha-profasi-definition.png) |
|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 6**: The instruction (14) -pappu- creates the left image and the (15) -profasi- creates the rigth image.

| ![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-custom-region-customized-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-custom-region-customized-definition.png)|
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 7**: Both instructions (16) and (17) create this image.

Note:
If the user, when using the **customized\_region** or **cr** option, specifies another conformation the program will generate two images: one figure will correspond to the custom region and conformation, and the other will be blank.
```
$ ./ss-map.py the/path/to/the/numpy_array.npy -a -cr conf 0 30 -170 120                 (18)
```
| ![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-both-estructura-alpha-customized-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-both-estructura-alpha-customized-definition.png) | ![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-both-estructura-custom-region-customized-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-both-estructura-custom-region-customized-definition.png)|
|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 8**: Both images are generated at the same time by the instruction (18). The left image is the one generated due to the -a argument, the right image is the one generated by the -cr argument.


## Calling the STRIDE ##

The program can call the program STRIDE independently from the other options.
```
$ ./ss-map.py the/path/to/the/numpy_array.npy -stride /path/to/the/pdb/files alpha       (19)
or 
$ ./ss-map.py the/path/to/the/numpy_array.npy -st /path/to/the/pdb/files beta            (20)
```
| ![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-stride-alpha.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-estructura-stride-alpha.png)| ![https://ss-map.googlecode.com/svn/wiki/images/examples/beta-estructura-stride-beta.png](https://ss-map.googlecode.com/svn/wiki/images/examples/beta-estructura-stride-beta.png)|
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 9**: The instruction (19) generates the left image and the instruction (20) generates the right image.

Both conformation can be called at the same time and the instruction will generate two images at the same time, one for the alpha helix conformation and another for the beta-strand conformation.


## Changing the images ##

The following options allow the user to set the image properties.

By default the program generates images using the colormap jet defined by the library pylab, and the scale is automatically fitted from zero to the maximum percentage present in the data. The user can set the program to generate black and white images using the option **colormap**, and the minimum and maximum in the scale with the option **rang\_colorbar**.
```
$ ./ss-map.py the/path/to/the/numpy_array.npy -a -cm binary                              (21)

$ ./ss-map.py the/path/to/the/numpy_array.npy -a -rgc 0 0.5                              (22)

# -cm is the abreviature for -colormap and -rgc is the abreviature for -range_colorbar
```
|![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-bw-estructura-alpha-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-bw-estructura-alpha-blackledge-definition.png) | ![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-rgc-estructura-alpha-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-rgc-estructura-alpha-blackledge-definition.png)|
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 10**: The instruction (21) generates the left image and the instruction (22) generates the right image. If the options where the default ones both images would be identical to the left image in Figure 1.

Maybe just a region of all the figure has an interest to the user, in this case the user can specify the region to show in the image using the instructions **-groups** and **-residues**. The first instruction will set the y axe limits and the second one the x axe limits.

```
$ ./ss-map.py the/path/to/the/numpy_array.npy -a -g 28 39     (23)

$ ./ss-map.py the/path/to/the/numpy_array.npy -a -r 2 9       (24)
```
|![https://ss-map.googlecode.com/svn/wiki/images/examples/beta-g-estructura-alpha-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/beta-g-estructura-alpha-blackledge-definition.png) | ![https://ss-map.googlecode.com/svn/wiki/images/examples/beta-r-estructura-alpha-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/beta-r-estructura-alpha-blackledge-definition.png) |
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 11**: If these images were created using the default settings they would be identical to the rigth image in Figure 1.

## Other plots ##
The following instructions allow the user to generate other figures.

The **-helix-per-residue** option generates a figure representing the total alpha-helix percentage in the ensemble per residue.
```
$ ./ss-map.py the/path/to/the/numpy_array.npy -a -helix-per-residue  (25)
or
$ ./ss-map.py the/path/to/the/numpy_array.npy -a -hr                 (26)
```
|![https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-alpha-helix-percentage-per-residue-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/alpha-alpha-helix-percentage-per-residue-blackledge-definition.png) |
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 12**: In this figure the x axe represents the residue and the y axe represents the total percentage of alpha helix in the ensemble.

The **-helix-per-group** option represents the alpha-helix structure length at diferent temperatures (or ensembles).
```
$ ./ss-map.py the/path/to/the/numpy_array.npy  -helix-per-group      (27)
or 
$ ./ss-map.py the/path/to/the/numpy_array.npy  -hgt                  (28)
```
|![https://ss-map.googlecode.com/svn/wiki/images/examples/multi-helix-default-alpha-helix-percentage-per-group-and-temperature-blackledge.png](https://ss-map.googlecode.com/svn/wiki/images/examples/multi-helix-default-alpha-helix-percentage-per-group-and-temperature-blackledge.png) |
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 13**: In this figure the axe y represent the temperature and the x axe represents the structured region lenght.

The **-temperature** option allows the user to set the y axe values when using the **-hgt** option.
```
$ ./ss-map.py the/path/to/the/numpy_array.npy  -hgt -temperature 335 327 320 313 306 299 292 286 280 274     (29)
or
$ ./ss-map.py the/path/to/the/numpy_array.npy  -hgt -temp 335 327 320 313 306 299 292 286 280 274            (30)
```
|![https://ss-map.googlecode.com/svn/wiki/images/examples/multi-helix-temperatures-alpha-helix-percentage-per-group-and-temperature-blackledge.png](https://ss-map.googlecode.com/svn/wiki/images/examples/multi-helix-temperatures-alpha-helix-percentage-per-group-and-temperature-blackledge.png)|
|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 14**: This image is the same as the one in **Figure 13** but specifying the y axes values.

## Saving the information ##
The following instructions allow the user to save the data generated by the program.

The option **-save\_figure** or **-sf** allows the user to save automatically the figure generated. The figure is saved with the following structure: **figure\_name-**_alpha-helix/beta-strand/poliprolineII-helix/custom-region_**-structure-**_blackledge/pappu/profasi/customized_**-definition.png**
```
$ ./ss-map.py the/path/to/the/numpy_array -a -sf path/to/save/the/figure/figure_name    (29)
```

The option **-save\_numpy** or **-sn** allows the user to save the percentages in a numpy array. The numpy is saved using the following structure: numpy\_name-_alpha-helix/beta-strand/poliprolineII-helix/custom-region_-percentage-_blackledge/pappu/profasi/customized_-definition.npy
```
$ ./ss-map.py the/path/to/the/numpy_array -a -sf path/to/save/the/figure/numpy_name    (30)
```

The option **-txt** will create multiple files, precisely it will create one .txt file for each structure studied. This files will contain the one percentage array for file. The files names will have the structure: **ss-map-**_(stride)/alpha/beta/poliprolyne/custom-region_**-percentage-**_blackledge/pappu/profasi/customized_**-definition.txt**

## Combining Options ##

All the options can be combined, this will create as many images as conformations are selected.

```
$ ./ss-map.py the/path/to/the/numpy_array -a -b -st path/to/the/pdbs/ beta -sf /path/to/save/the/figure/general_name   (31)
```
|![https://ss-map.googlecode.com/svn/wiki/images/examples/combining-options-alpha-helix-estructure-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/combining-options-alpha-helix-estructure-blackledge-definition.png)|![https://ss-map.googlecode.com/svn/wiki/images/examples/combining-options-beta-strand-estructure-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/combining-options-beta-strand-estructure-blackledge-definition.png)|![https://ss-map.googlecode.com/svn/wiki/images/examples/combining-options-stride-beta-strand-estructure-blackledge-definition.png](https://ss-map.googlecode.com/svn/wiki/images/examples/combining-options-stride-beta-strand-estructure-blackledge-definition.png)|
|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
**Figure 15**: The three images are created at the same time using the instruction (31).

The three images are saved in the same directory, _/path/to/save/the/figure/_, and the names will be **_general\_name_-alpha-helix-estructure-blackledge-definition.png**, **_general\_name_-beta-strand-estructure-blackledge-definition.png and _general\_name_-stride-beta-strand-estructure-blackledge-definition.png**. In this instruction the structure\_definition option isn't used so the default option (blackledge) is being used.

