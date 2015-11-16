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

files = sys.argv[1:]
traj = md.load(files)
phi = md.geometry.compute_phi(traj)[1] #from -pi to pi
psi = md.geometry.compute_psi(traj)[1] #from -pi to pi
angles = np.rollaxis(np.array([psi, phi]), 0, start=3)
angles *=180/np.pi

np.save('angles', angles)
