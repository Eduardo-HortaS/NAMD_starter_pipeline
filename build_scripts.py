#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  build_scripts.py
#  
#  Copyright 2020 Fernando <fernando@winter>
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




from NaMD_Object3 import *
from pdb_to_pdb_fixed import export_PDB_fixed

export_PDB_fixed()

crd           = 'complex.inpcrd'  
parm          = 'complex.prmtop'
pdb_fixed1    = 'complex_fix.pdb'



outputname_00 = 'complex_solvent_equilibration'
outputname_01 = 'complex_system_equilibration'
outputname_02 = 'complex_production'



# creating NaMD object 
namdObject = NamdObject(sysType        = 'amber'      ,
                        coords         = crd          ,
                        parm           = parm         ,
                        psf            = None         ,
                        fixedAtomsFile = pdb_fixed1   )





# Solvent Equilibration 
namdObject.fixedAtoms   = True
namdObject.ZeroMomentum = False
namdObject.build_NaMD_input2(
	                         outputname = outputname_00,
	                         restart    = False        ,
	                         minimize   = 1000         ,
	                         dynamics   = 10000        ,
	                         parm       =  None        , 
	                         ensemble   = 'npt' )


# Full sistem Equilibration 
namdObject.fixedAtoms   = False
namdObject.binCoordinates = outputname_00+'.coor'
namdObject.build_NaMD_input2(
	                         outputname = outputname_01,
	                         restart    = False        ,
	                         minimize   = 1000         ,
	                         dynamics   = 10000        ,
	                         parm       =  None        , 
	                         ensemble   = 'npt' )


# Data collection
namdObject.fixedAtoms   = False
namdObject.binCoordinates = outputname_01+'.coor'
#namdObject.binVelocities  = outputname_01+'.vel'
#namdObject.extendedSystem = outputname_01+'.xsc'
namdObject.build_NaMD_input2(
	                         outputname = outputname_02,
	                         restart    = False        ,
	                         minimize   = 1000         ,
	                         dynamics   = 10000        ,
	                         parm       =  None        , 
	                         ensemble   = 'npt' )

