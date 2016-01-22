#!/usr/bin/env python3
#
#============================ SS-map ==================================
# get_angles is a Python program to generate a list of phi-psi angles
# from different structure files (using mdtraj).
# It is used to generate the input of ss-map
# get_angles is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# get_angles is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with ss-map.  If not, see <http://www.gnu.org/licenses/>.
# Written by: Jelisa Iglesias.
#
#======================================================================
#
# Author:
# Ramon Crehuet
#======================================================================
# For more information visit:
# https://github.com/rcrehuet/ss-map
# This work is done by the Theoretical and Computational Group
# of the IQAC (CSIC)
# http://iqac.csic.es/qteor
#
import sys
import argparse

try: 
    import numpy as np
except ImportError:
    print("get_angle needs numpy library.")
    sys.exit()

try:
    import mdtraj as md
except ImportError:
    print("get_angles needs the mdtraj library. Get it from http://mdtraj.org")
    sys.exit()


parser = argparse.ArgumentParser(\
      description="Generate an numpy array with the phi and psi angles.")
parser.add_argument("filename", nargs='+', 
    help="The structure file(s) (pdb, xtc,...).")
parser.add_argument('--top', 
                   help='The topology file (needed for structure files \
                   without topology)')
args = parser.parse_args()

if not args.top:
    try:
        traj = md.load(args.filename)
    except ValueError:
        print("ERROR: You need to specify a topology (pdb, gro, ...) for this \
kind of trajectory file.")
        sys.exit()
else:
    traj = md.load(args.filename, top=args.top)

atoms_phi, phi = md.geometry.compute_phi(traj) #from -pi to pi
atoms_psi, psi = md.geometry.compute_psi(traj) #from -pi to pi
atoms_phi = atoms_phi[:, 1]
atoms_psi = atoms_psi[:, 0]
phi = phi[:, np.in1d(atoms_phi, atoms_psi)]
psi = psi[:, np.in1d(atoms_psi, atoms_phi)]

angles = np.rollaxis(np.array([psi, phi]), 0, start=3)
angles *=180/np.pi

np.save('angles', angles)
