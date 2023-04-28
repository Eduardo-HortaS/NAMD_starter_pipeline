#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pdb_to_pdb_fixed.py
#  
#  Copyright 2020 Fernando <fernando@Frost>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import os, sys

def export_PDB_fixed (pdbin   = 2h1x_restrict.pdb, 
                      pdbout  = fix.pdb, 
                      wat     = True,
                      lig     = False,
                      protein = False, 
                      
                      ):
	""" Function doc """
	data = open(pdbin, 'r')
	data  =  data.readlines()
	
	data2 = []
	for line in data:
		
		if 'WAT' in line:
			line = line.replace('  1.00  ','  0.00  ')
		
		data2.append(line)
		
	fileout = open(pdbout, 'w')
	fileout.writelines(data2)


print sys.argv
args = sys.argv

pdbin = args[1]

pdbout = args[2]

export_PDB_fixed(pdbin = pdbin, pdbout = pdbout)
