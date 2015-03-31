# Frequently Asked Questions #

## Why the first and the last aminoacids in a pdb file aren't studied? ##

> The first aminoacid does not present a phi angle and the last aminoacid does not present a psi angle, due to the lack of a previous and a next aminoacid respectively. Thus they cannot be situated in the Ramachandran space.

## Errors ##
  * Why does the instruction _./ss-map.py  -a -st path/to/the/pdbs/ beta the/path/to/the/numpy\_array_ give an error?
    * It is because the option **-stride** or **-st** takes all the following arguments, thus taking _the/path/to/the/numpy\_array_ as part of the stride argument and not as the only mandatory parameter, the path to the array.

