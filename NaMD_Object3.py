#import numpy as np
from pprint import pprint
import os
import sys
from pprint import pprint

print(sys.argv)
args = sys.argv

#filein = open('/home/fernando/Documents/Betao/sug1.BL03000001_box.pdb', 'r')



class NamdObject:
	""" Class doc """
	
	def __init__ (self, sysType         = 'amber', 
	                     coords         = None   , 
	                     parm           = None   , 
	                     psf            = None   , 
	                     fixedAtomsFile = None   ):
		
		
		""" Class initialiser """
		
		# NamdObject atributes 
		
		self.minimize         = True                               
		self.interactions     = 1000                               
		self.dynamics         = False                              
		self.nStep            = 10000                              
		
		if sysType == 'amber':
			self.amber            = True # yes / no
			self.parmfile         = parm #topin  
		else:
			self.amber            = False
			self.parmfile         = None #topin  

		if sysType == 'charmm':
			self.paraTypeCharmm	  = True
			self.parameters       = parm #  'par_all27_prot_lipid.inp'
			self.psf              = psf #  structure mypsf.psf
		                      
		self.coordinates      = coords #pdbin                              
		
		#self.coord_type       = 'crd'                              

		self.firsttimestep    = 0
		
		#If we are starting from scratch, we'll use the coordinates from the PDB file and take velocities randomly from a Boltzmann distribution, using the $temperature variable.
		# starting from scratc
		self.temperature      = 300                                
		
		
		#Otherwise we'll read in the binary output files of the previous run. We use the Tcl variable $inputname to avoid errors in typing the input file names, since we will end up copying and modifying this file several times during the course of a long simulation. The extendedSystem file holds the periodic cell dimensions that are needed to continue after a constant pressure simulation. In order to start numbering timesteps with the final step of the previous run, we use firsttimestep to specify the timestep on which the input coordinates and velocities were written. The number of steps calculated in this run will be numsteps - firsttimestep.
		# continuing a run
		self.binCoordinates   = None              #$inputname.restart.coor                                                        | coordinates from last run (binary) 
		self.binVelocities    = None              #$inputname.restart.vel  ;# remove the "temperature" entry if you use this!     | velocities from last run (binary) 
		self.extendedSystem	  = None              #$inputname.xsc                                                                 | cell dimensions from last run 
		self.firsttimestep    = 0                 # last step of previous run                                                     | last step of previous run 
		self.numsteps         = 100000            # run stops when this step is reached                                           | run stops when this step is reached 
		
		#When initially assembling a system it is sometimes useful to fix the protein while equilibrating water or lipids around it. These options read a PDB file containing flags for the atoms to fix. The number and order of atoms in the PDB file must exactly match that of the PSF and other input files.
		self.fixedAtoms       = False                              
		self.fixedAtomsFile   = fixedAtomsFile  
		self.fixedAtomsCol    = None #  A/B                 ;# set beta non-zero to fix an atom

		
		self.readexclusions  = True
		self.scnb            = 2.0

		# Force-Field Parameters
		# These are specified by CHARMM
		self.exclude        = 'scaled1-4'
		self.scaling1_4     = 1.0
		self.cutoff         = 12.0 # may use smaller, maybe 10., with PME
		self.switching      = True 
		self.switchdist     = 10.0 # cutoff - 2.
		self.pairlistdist   = 14.0 # cutoff + 2.

		
		# These are specified by AMBER
		self.pairListDist       =	11            
		self.LJcorrection       =	True          
		self.ZeroMomentum       =	True          
		self.rigidBonds         =	'{water, all}'
		self.rigidTolerance     =	1.0e-8        
		self.rigidIterations    =	100           
		self.useSettle          =	True          
		self.PMETolerance       =	1.0e-6
		self.PMEInterpOrder     =	4     


		# Integrator Parameters
		self.timestep           = 2.0    # 2fs/step
		self.rigidBonds         = 'all'  # needed for 2fs steps
		self.nonbondedFreq      = 1      # nonbonded forces every step
		self.fullElectFrequency = 2      #
		
		#The integration timestep is normally limited to 1fs. This can be extended to 2fs by fixing the length of all bonds in the molecule involving hydrogen atoms via rigidBonds all, which also makes water molecules completely rigid. To only make water rigid, use rigidBonds water and a 1fs timestep. For any simulations involving water molecules, one should make water rigid since water molecules have been parametrized as rigid water. Nonbonded forces must be calculated at least every 2fs, which in this example is every step. Full electrostatics forces (from particle mesh Ewald, discussed below) are evaluated every other step in this example. nonbondedFreq and fullElectFrequency must evenly divide stepspercycle.
		# Integrator Parameters
		self.timestep           = 2.0    # 2fs/step

		#The chosen stepspercycle is related to pairlistdist above, as pairlists are generated at the beginning of each cycle and are assumed to contain every atom that passes within the cutoff distance during entire cycle. Warnings will appear in the output if pairlistdist is too short (or stepspercycle is too large). The length of the simulation must be an integer number of cycles.
		self.stepspercycle      = 10     # redo pairlists every ten steps


		# Langevin dynamics balances friction with random noise to drive each atom in the system towards a target temperature. The following parameters are good for equilibration, but a langevinDamping as low as 1. would be sufficient for simulation. Higher damping may be necessary if, for example, work done on the system by steering forces is driving up the temperature.
		# Constant Temperature Control
		self.langevin           = True #on      # langevin dynamics
		self.langevinDamping    = 1.0           # damping coefficient of 1/ps
		#self.langevinTemp       = $temperature # random noise at this level
		self.langevinHydrogen   = False         # don't couple bath to hydrogens

		

		#We will write restart files, coordinate trajectory .dcd files, and extended system (periodic cell) trajectory .xst files at regular intervals.
		self.restartfreq      = 1000                               
		self.dcdfreq          = 1000                               
		self.xstFreq          = 1000
		self.outputEnergies   = 1000                               
		self.outputPressure   = 1000                               						
		self.outputTiming     = 1000     ;# shows time per step and time to completion
		
		#The outputName prefix will be used to create all of the trajectory (.dcd, .xst), output (.coor, vel, .xsc), and restart (.restart.coor, .restart.vel, .restart.xsc) files generated by this run. Be careful that this is different from $myinput, or the job will overwrite its input files!
		self.outputname       = 'system_minimized'                 
		
		
		self.cellBasisVector1 = None  #cell_parameters['cellBasisVector1']
		self.cellBasisVector2 = None  #cell_parameters['cellBasisVector2']
		self.cellBasisVector3 = None  #cell_parameters['cellBasisVector3']
		self.cellOrigin       = None  #cell_parameters['cellOrigin'      ]
		self.wrapAll          = True                               
		
		self.wrapWater        = True  #on              ;# wrap water to central cell
		self.wrapAll          = True  #on              ;# wrap other molecules too
		self.wrapNearest      = False #off             ;# use for non-rectangular cells
 	


		#Particle mesh Ewald (PME) full electrostatics are more accurate and less expensive than larger cutoffs, and are recommended for most work. PME is only applicable to periodic simulations, and the user must specify a grid corresponding to the size of the periodic cell. PME grid dimensions should have small integer factors only and be greater than or equal to length of the basis vector. To manually define the grid size instead of letting NAMD choose the dimensions for you according to these guidelines, replace the GridSpacing command with explicit GridSize commands instead.
		#PME (for full-system periodic electrostatics)
		self.PME              = True 
		self.PMEGridSpacing   = 1.0                                     
		self.PMEGridSizeX     = None #cell_parameters['PMEGridSizeX']  # 2^5, close to 31.2               
		self.PMEGridSizeY     = None #cell_parameters['PMEGridSizeY']  # 3^2 * 5, close to 44.8              
		self.PMEGridSizeZ     = None #cell_parameters['PMEGridSizeZ']  # 2 * 3^3, close to 51.3              

		#Constant pressure is recommended for periodic simulations. Using group-based pressure to control the periodic cell fluctuations is desirable because the atom-based pressure has more high-frequency noise. useFlexibleCell is useful for anisotropic systems such as membranes, allowing the height, length, and width of the cell to vary independently, possibly even fixing the lipid cross-sectional (x-y plane) area with useConstantArea. For a protein surrounded by water there is nothing to prevent the cell from becoming highly extended in one dimension, so it is better to choose useFlexibleCell no in this case.
		# Constant Pressure Control (variable volume)
		self.useGroupPressure = True    # yes ;# needed for rigid bonds
		self.useFlexibleCell  = False   # no  ;# no for water box, yes for membrane
		self.useConstantArea  = False   # no  ;# no for water box, maybe for membrane


		#The actual parameters of the Nose-Hoover Langevin piston method control the target pressure and the dynamical properties of the barostat. A ``piston'' with a longer period (i.e., larger mass) will better damp out fluctuations in the instantaneous pressure. Langevin dynamics is applied to the piston itself coupling it to a heat bath with a damping constant of $1 / $langevinPistonDecay. We set langevinPistonDecay smaller that langevinPistonPeriod to ensure that harmonic oscillations in the periodic cell are overdamped.
		self.langevinPiston        = True #on
		self.langevinPistonTarget  = 1.01325      ;# pressure in bar -> 1 atm
		self.langevinPistonPeriod  = 100.         ;# oscillation period around 100 fs
		self.langevinPistonDecay   = 50.          ;# oscillation decay time of 50 fs
		#self.langevinPistonTemp    $temperature ;# coupled to heat bath

		
		##The interactive MD features of NAMD and VMD allow you to connect to a running simulation to apply steering forces manually. These options affect performance, and should therefore not be used unless you are actually steering the simulation.
		#self.IMDon        =   True       #on
		#self.IMDport      =   3000       # port number (enter it in VMD)
		#self.IMDfreq      =   1          # send every 1 frame
		#self.IMDwait      =   True #no   # wait for VMD to connect before running?
		
		#Now we minimize the system to eliminate bad initial contacts, reinitialize the velocities to the desired target temperature (since minimization sets velocities to zero), and run for 100ps. We could accomplish the same thing with two different NAMD runs using the numsteps and minimization options. Scripting commands such as those below override numsteps.
		#self.minimize            1000          ;# lower potential energy for 1000 steps
		#self.reinitvels          $temperature  ;# since minimization zeros velocities
		#self.run 50000 ;# 100ps

		self.get_cell ()
		if sysType == 'amber':
			self.set_amberparm ()
		
		self.namd_path  = '/home/peagah/NAMD_2.14_Linux-x86_64-multicore/namd2'
		self.nCPUs      = 12
		self.nGPUs      = 0


	def set_amberparm (self, parameters = {} ):
		""" Function doc """
		self.amber              =   True                
		self.switching          =	False           #off           
		self.exclude            =	'scaled1-4'     #scaled1-4     
		self.scaling1_4        	=   0.833333333     #0.833333333   
		self.scnb               =	2.0             #2.0           
		self.readexclusions     =	True            #yes           
		self.cutoff             =	9               #9             
		self.watermodel         =	'tip3'          #{tip3,        
		self.pairListDist       =	11              #11            
		self.LJcorrection       =	True            #on            
		self.ZeroMomentum       =	True            #on            
		self.rigidBonds         =	'all'           # {water, all}     
		self.rigidTolerance     =	1.0e-8          #1.0e-8        
		self.rigidIterations    =	100             #100           
		self.useSettle          =	True            #{on}          
		self.timeStep           =	1.0             #1.0           
		
		
		# Non-bonded interactions 
		self.fullElectFrequency =	1               #1           
		self.nonBondedFreq      =	1               #1             
		self.stepspercycle      =	10              #10            
		#self.PME               =	True            #on            
		
		#self.PMEGridSizeX       =	intx,           #intx,
		#self.PMEGridSizeY       =	inty,           #inty,
		#self.PMEGridSizeZ       =	intz            #intz          
		self.PMETolerance       =	1.0e-6          #1.0e-6        
		self.PMEInterpOrder     =	4               #4             
		#self.cellBasisVector1   =	x1              #x1            
		#self.cellBasisVector1   =	x2              #x2            
		#self.cellBasisVector1   =	x3              #x3            
		#self.cellOrigin         =	x               #x             
		
		
		self.langevin                        =  True       # {on,off} True/False         
		self.langevinDamping                 =  5          # 5          
		self.langevinTemp                    =  300        # 300        
		self.langevinHydrogen                =  False      # off        
 
		self.BerendsenPressure               =  False   # {on,off} True/False 
		self.BerendsenPressureTarget         =  1.01325    # 1.01325    
		self.BerendsenPressureRelaxationTime =  100.0      # 100.0      
		
		self.useGroupPressure                =  True        # yes
		self.useFlexibleCell                 =  False         # no
		self.useConstantArea                 =  False         # no         
		

	def  build_NaMD_input2( self                  , 
	                        outputname = 'test'   ,
	                        restart    = False    ,
	                        minimize   = 1000     ,
	                        dynamics   = 1000000  ,
	                        parm       =  None    , 
	                        ensemble   = 'npt'    ):
	
		restart = False
	
		configFile = open(outputname + '.namd', 'w')
		
		text = ''
	
		text += '#############################################################                \n'
		text += '##                   JOB DESCRIPTION                       ##                \n'
		text += '#############################################################              \n\n'
		
		#text += '# Minimization and                                                         \n\n'
		#
		#text += '#############################################################                \n'
		#text += '## ADJUSTABLE PARAMETERS                                   ##                \n'
		#text += '#############################################################              \n\n'
		
		#text += '#reiniciando a simulacao                                                     \n'
		#text += '#bincoordinates     .restart.coor                                            \n'
		#text += '#binvelocities      .restart.vel                                             \n'
		#text += '#extendedSystem     .restart.xsc                                             \n'
		#text += '#structure          .psf                                                     \n'

		if self.amber == True:
			text += '\n\n'
			text += '#############################################################                \n'
			text += '##                   amber simulation                      ##                \n'
			text += '#############################################################                \n\n'
			text += 'amber                           yes    \n'
			text += 'parmfile                        {:100} \n'  .format(self.parmfile)
		
		else:
			if paraTypeCharmm == True:
				
				text += '\n\n'
				text += '#############################################################                \n'
				text += '##                 charmm  simulation                      ##                \n'
				text += '#############################################################                \n'
				text += 'paraTypeCharmm                  yes \n '
			else:
				text += 'paraTypeCharmm                  no \n '
			text += 'parameters                      {:100} \n'  .format(self.parameters)
			text += 'structure                       {:100} \n'  .format(self.psf)
		

		
		#defining coordinate type
		coordType = self.coordinates.split('.')
		
		if coordType[-1] == 'crd' or coordType[-1] == 'inpcrd' or coordType[-1] == 'rst':
			
			text += 'ambercoor                       {:<100} \n' .format(self.coordinates) 

		if coordType[-1] == 'pdb':
			text += 'coordinates                     {:<100} \n' .format(self.coordinates)
		

		if self.binCoordinates == None:
			pass
		else:               
			text += 'bincoordinates                  {:<100} \n' .format(self.binCoordinates)	

		
		#  fixedAtoms ?
		if self.fixedAtoms == True:
			text += '\n\n# mobile atom selection:   \n'
			text += 'fixedAtoms                      on     \n'
			text += 'fixedAtomsFile                  {:<100} \n' .format(self.fixedAtomsFile)	
			
			if self.fixedAtomsCol == None:
				pass
			else:
				text += 'fixedAtomsCol                   {:<100} \n' .format(self.fixedAtomsCol)

		text += 'set outputname                  {:<100} \n' .format(outputname)
		
		#---------------------------------------------------------------------------------------------------------------------------------------#
		#              - - - - - - - - - - - - - - - - - - - -   R E S T A R T   - - - - - - - - - - - - - - - - - - - - - - -                  #
		#---------------------------------------------------------------------------------------------------------------------------------------#
		if self.binVelocities != None:#self.restart == True:                                                                                    #
			if self.extendedSystem != None:                                                                                                     #
				text += 'binVelocities                   {:<100} \n'.format (self.binVelocities )                                               #
				text += 'extendedSystem                  {:<100} \n'.format (self.extendedSystem)                                               #
				text += 'firsttimestep                   {:<100} \n'.format (self.firsttimestep)                                                #
																																				#
				if self.readexclusions == True:                                                                                                 #
					text += 'readexclusions                  yes    \n' #self.scnb                                                              #
					text += 'scnb                            {:<5}  \n\n'.format(self.scnb)		                                                #
																																				#
				#''' - - - - - - - - -     A M B E R  T Y P E      - - - - - - - - - - -'''                                                      #
				if self.amber == True:                                                                                                          #
				#text += 'amber                           {:<5} \n'.format(self.amber          ) #   =   True                                   #
					text += '#############################################################                \n'                                   #
					text += '## SIMULATION PARAMETERS                                   ##                \n'                                   #
					text += '#############################################################              \n\n'                                   #
					if self.switching == False:                                                                                                 #
						text += 'switching                       off \n'                                                                        #
					else:                                                                                                                       #
						text += 'switching                       on \n'                                                                         #
						text += 'switchdist                      {:<20} \n'.format(self.switchdist  ) #    10.                                  #
						text += 'pairlistdist                    {:<20} \n'.format(self.pairlistdist) #    13.5                                 #
																																				#
																																				#
					text += 'exclude                         {:<20} \n'.format(self.exclude        ) #   =	'scaled1-4'     #scaled1-4          #
					text += '1-4scaling                      {:<20} \n'.format(self.scaling1_4     ) #  	=   0.833333333     #0.833333333    #
					text += 'cutoff                          {:<20} \n'.format(self.cutoff         ) #   =	9               #9                  #
					text += 'watermodel                      {:<20} \n'.format(self.watermodel     ) #   =	'tip3'          #{tip3,             #
					text += 'pairListDist                    {:<20} \n'.format(self.pairListDist   ) #   =	11              #11                 #
					if self.LJcorrection == True:                                                                                               #
						text += 'LJcorrection                    on \n'                                                                         #
					if self.ZeroMomentum == True:                                                                                               #
						text += 'ZeroMomentum                    on \n'                                                                         #
																																				#
					text += 'rigidBonds                      {:<20} \n'.format(self.rigidBonds     ) #   =	'{water, all}'  # {water, all}      #
					text += 'rigidTolerance                  {:<20} \n'.format(self.rigidTolerance ) #   =	1.0e-8          #1.0e-8             #
					text += 'rigidIterations                 {:<20} \n'.format(self.rigidIterations) #   =	100             #100                #
																																				#
					if self.useSettle == True:                                                                                                  #
						text += 'useSettle                       on \n'.format(self.useSettle         )                                         #
					text += 'fullElectFrequency              {:<20} \n'.format(self.fullElectFrequency)                                         #
					text += 'nonBondedFreq                   {:<20} \n'.format(self.nonBondedFreq     )                                         #
					text += 'stepspercycle                   {:<20} \n'.format(self.stepspercycle     )                                         #
					text += 'timeStep                        {:<20} \n'.format(self.timeStep          )		                                    #
																																				#
																																				#
				#''' - - - - - - - - -     C H A R M M   T Y P E      - - - - - - - - - - -'''                                                   #
				else:                                                                                                                           #
					text += '#############################################################                \n'                                   #
					text += '## SIMULATION PARAMETERS                                   ##                \n'                                   #
					text += '#############################################################              \n\n'                                   #
					text += '# Input                                                                      \n'                                   #
					#       'amber                           yes    \n'                                                                         #
					text += '\n\n# Force-Field Parameters \n'                                                                                   #
					text += 'exclude                         {:<20} \n'.format(self.exclude     ) #    scaled1-4                                #
					text += '1-4scaling                      {:<20} \n'.format(self.scaling1_4  ) #    1.0                                      #
					text += 'cutoff                          {:<20} \n'.format(self.cutoff      ) #    12.                                      #
					text += 'switching                       {:<20} \n'.format(self.switching   ) #    on                                       #
					text += 'switchdist                      {:<20} \n'.format(self.switchdist  ) #    10.                                      #
					text += 'pairlistdist                    {:<20} \n'.format(self.pairlistdist) #    13.5                                     #
																																				#
					text += '\n\n# Integrator Parameters  \n'                                                                                   #
					text += 'timestep                        {:<20} \n'.format(self.timestep          )                                         #
					text += 'rigidBonds                      {:<20} \n'.format(self.rigidBonds        )                                         #
					text += 'nonbondedFreq                   {:<20} \n'.format(self.nonbondedFreq     )                                         #
					text += 'fullElectFrequency              {:<20} \n'.format(self.fullElectFrequency)                                         #
					text += 'stepspercycle                   {:<20} \n'.format(self.stepspercycle     )                                         #
																																				#
		#---------------------------------------------------------------------------------------------------------------------------------------#
		
		
		
		
		
		#---------------------------------------------------------------------------------------------------------------------------------------#
		#              - - - - - - - - - - - - - - - -   N E W   S I M U L A T I O N   - - - - - - - - - - - - - - - - - - - -                  #
		#---------------------------------------------------------------------------------------------------------------------------------------#
		else:                                                                                                                                   #
			text += 'set temperature                 {:<10}  \n' .format(self.temperature)                                                      #
			text += 'temperature                     $temperature \n'                                                                           #
			text += 'firsttimestep                   {:<100} \n' .format(self.firsttimestep)                                                    #
																																				#
			if self.readexclusions == True:                                                                                                     #
				text += 'readexclusions                  yes    \n' #self.scnb                                                                  #
				text += 'scnb                            {:<5}  \n\n'.format(self.scnb)                                                         #
																																				#
			#''' - - - - - - - - -     A M B E R  T Y P E      - - - - - - - - - - -'''                                                         #
			if self.amber == True:                                                                                                              #
			#text += 'amber                           {:<5} \n'.format(self.amber          ) #   =   True                                       #
				text += '#############################################################                \n'                                       #
				text += '## SIMULATION PARAMETERS                                   ##                \n'                                       #
				text += '#############################################################              \n\n'                                       #
				if self.switching == False:                                                                                                     #
					text += 'switching                       off \n'                                                                            #
				else:                                                                                                                           #
					text += 'switching                       on \n'                                                                             #
					text += 'switchdist                      {:<20} \n'.format(self.switchdist  ) #    10.                                      #
					text += 'pairlistdist                    {:<20} \n'.format(self.pairlistdist) #    13.5                                     #
																																				#
																																				#
				text += 'exclude                         {:<20} \n'.format(self.exclude        ) #   =	'scaled1-4'     #scaled1-4              #
				text += '1-4scaling                      {:<20} \n'.format(self.scaling1_4     ) #  	=   0.833333333     #0.833333333        #
				text += 'cutoff                          {:<20} \n'.format(self.cutoff         ) #   =	9               #9                      #
				text += 'watermodel                      {:<20} \n'.format(self.watermodel     ) #   =	'tip3'          #{tip3,                 #
				text += 'pairListDist                    {:<20} \n'.format(self.pairListDist   ) #   =	11              #11                     #
				if self.LJcorrection == True:                                                                                                   #
					text += 'LJcorrection                    on \n'                                                                             #
				if self.ZeroMomentum == True:                                                                                                   #
					text += 'ZeroMomentum                    on \n'                                                                             #
																																				#
				text += 'rigidBonds                      {:<20} \n'.format(self.rigidBonds     ) #   =	'{water, all}'  # {water, all}          #
				text += 'rigidTolerance                  {:<20} \n'.format(self.rigidTolerance ) #   =	1.0e-8          #1.0e-8                 #
				text += 'rigidIterations                 {:<20} \n'.format(self.rigidIterations) #   =	100             #100                    #
																																				#
				if self.useSettle == True:                                                                                                      #
					text += 'useSettle                       on \n'.format(self.useSettle         )                                             #
				text += 'fullElectFrequency              {:<20} \n'.format(self.fullElectFrequency)                                             #
				text += 'nonBondedFreq                   {:<20} \n'.format(self.nonBondedFreq     )                                             #
				text += 'stepspercycle                   {:<20} \n'.format(self.stepspercycle     )                                             #
				text += 'timeStep                        {:<20} \n'.format(self.timeStep          )		                                        #
																																				#
																																				#
			#''' - - - - - - - - -     C H A R M M   T Y P E      - - - - - - - - - - -'''                                                      #
			else:                                                                                                                               #
				text += '#############################################################                \n'                                       #
				text += '## SIMULATION PARAMETERS                                   ##                \n'                                       #
				text += '#############################################################              \n\n'                                       #
				text += '# Input                                                                      \n'                                       #
				#       'amber                           yes    \n'                                                                             #
				text += '\n\n# Force-Field Parameters \n'                                                                                       #
				text += 'exclude                         {:<20} \n'.format(self.exclude     ) #    scaled1-4                                    #
				text += '1-4scaling                      {:<20} \n'.format(self.scaling1_4  ) #    1.0                                          #
				text += 'cutoff                          {:<20} \n'.format(self.cutoff      ) #    12.                                          #
				text += 'switching                       {:<20} \n'.format(self.switching   ) #    on                                           #
				text += 'switchdist                      {:<20} \n'.format(self.switchdist  ) #    10.                                          #
				text += 'pairlistdist                    {:<20} \n'.format(self.pairlistdist) #    13.5                                         #
																																				#
				text += '\n\n# Integrator Parameters  \n'                                                                                       #
				text += 'timestep                        {:<20} \n'.format(self.timestep          )                                             #
				text += 'rigidBonds                      {:<20} \n'.format(self.rigidBonds        )                                             #
				text += 'nonbondedFreq                   {:<20} \n'.format(self.nonbondedFreq     )                                             #
				text += 'fullElectFrequency              {:<20} \n'.format(self.fullElectFrequency)                                             #
				text += 'stepspercycle                   {:<20} \n'.format(self.stepspercycle     )                                             #
																																				#
																																				#
			if self.langevin == True:                                                                                                           #
				text += '\n\n# Constant Temperature Control                                               \n'                                   #
				text += 'langevin                        yes    ;# do langevin dynamics                            \n'                          #
				text += 'langevinDamping                 {:<5}   ;# damping coefficient (gamma) of 2/ps\n'.format(self.langevinDamping)         #
				text += 'langevinTemp                    $temperature                                             \n'                           #
				text += 'langevinHydrogen                off    ;# dont couple langevin bath to hydrogens         \n'                           #
																																				#
			#if self.BerendsenPressure  ==  False:          # {on,off} True/False                                                               #
			#	text += 'BerendsenPressureTarget         {:<20)'.format(self.BerendsenPressureTarget        )#    =  1.01325    # 1.01325       #
			#	text += 'BerendsenPressureRelaxationTime {:<20)'.format(self.BerendsenPressureRelaxationTime)#    =  100.0      # 100.0         #
																																				#
																																				#
			if ensemble == 'npt':                                                                                                               #
																																				#
				text += '\n\n# Constant Pressure Control (variable volume)                                \n'                                   #
				#           'amber                           yes    \n'                                                                         #
				if self.useGroupPressure == True:                                                                                               #
					text += 'useGroupPressure                yes \n'                                                                            #
				else:                                                                                                                           #
					text += 'useGroupPressure                no \n'                                                                             #
																																				#
				if self.useFlexibleCell  == False:                                                                                              #
					text += 'useFlexibleCell                 no  \n'                                                                            #
				else:                                                                                                                           #
					text += 'useFlexibleCell                 yes \n'                                                                            #
																																				#
				if  self.useConstantArea  == False:   # no  ;# no for water box, maybe for membrane                                             #
					text += 'useConstantArea                 no  \n'                                                                            #
				else:                                                                                                                           #
					text += 'useConstantArea                 yes  \n'                                                                           #
																																				#
				if self.langevinPiston    == True:                                                                                              #
					#       'amber                           yes    \n'                                                                         #
					text += 'langevinPiston                  yes  \n'                                                                           #
					text += 'langevinPistonTarget            {:<20}  \n'.format(self.langevinPistonTarget) #1.01325 ;#  in bar -> 1 atm         #
					text += 'langevinPistonPeriod            {:<20}  \n'.format(self.langevinPistonPeriod) #100.                                #
					text += 'langevinPistonDecay             {:<20}  \n'.format(self.langevinPistonDecay ) #50.                                 #
					text += 'langevinPistonTemp              $temperature \n'                                                                   #
																																				#
		#---------------------------------------------------------------------------------------------------------------------------------------#

		
		#---------------------------------------------------------------------------------------------------------------------------------------#
		#              - - - - - - - - - - - - - - - -   C E L L   P A R A M E T E R S   - - - - - - - - - - - - - - - - - - -                  #
		#---------------------------------------------------------------------------------------------------------------------------------------#
	
		text += '\n\n# Periodic Boundary Conditions                                               \n'
		text += 'cellBasisVector1    {:>8.5f}    {:>8.5f}    {:>8.5f}   \n'.format(self.cellBasisVector1[0], self.cellBasisVector1[1], self.cellBasisVector1[2])
		text += 'cellBasisVector2    {:>8.5f}    {:>8.5f}    {:>8.5f}   \n'.format(self.cellBasisVector2[0], self.cellBasisVector2[1], self.cellBasisVector2[2])
		text += 'cellBasisVector3    {:>8.5f}    {:>8.5f}    {:>8.5f}   \n'.format(self.cellBasisVector3[0], self.cellBasisVector3[1], self.cellBasisVector3[2])
		text += 'cellOrigin          {:>8.5f}    {:>8.5f}    {:>8.5f}   \n'.format(self.cellOrigin      [0], self.cellOrigin      [1], self.cellOrigin      [2])
		
		if self.wrapAll == True:
			text += 'wrapAll                         on     \n'
		else:                                               
			text += 'wrapAll                         off    \n'

		#---------------------------------------------------------------------------------------------------------------------------------------#
		#              - - - - - - - - - - - - - - - - -   P M E   P A R A M E T E R S   - - - - - - - - - - - - - - - - - - -                  #
		#---------------------------------------------------------------------------------------------------------------------------------------#
	
		if self.PME  == True:
			text += '\n\n# PME (for full-system periodic electrostatics)                              \n'
			text += 'PME                             yes                                                      \n'
			text += 'PMEGridSpacing                  {:<10}    \n'.format(self.PMEGridSpacing)
			text += 'PMEGridSizeX                    {:<10}    \n'.format(self.PMEGridSizeX  )
			text += 'PMEGridSizeY                    {:<10}    \n'.format(self.PMEGridSizeY  )
			text += 'PMEGridSizeZ                    {:<10}    \n'.format(self.PMEGridSizeZ  )
			text += 'PMETolerance                    {:<10}    \n'.format(self.PMETolerance  )     
			text += 'PMEInterpOrder                  {:<10}    \n'.format(self.PMEInterpOrder)



		#---------------------------------------------------------------------------------------------------------------------------------------#
		#              - - - - - - - - - - - - - - - - -  O U T P U T S   - - - - - - - - - - - - - - - - - - -                                 #
		#---------------------------------------------------------------------------------------------------------------------------------------#
	
		text += '\n\n# Output                                                 \n'
		text += 'outputName                      $outputname                          \n'
		text += 'restartfreq                     {:<10}     ;# 1000steps = every 2 ps  \n'   .format(self.restartfreq   )
		text += 'dcdfreq                         {:<10}                                \n'   .format(self.dcdfreq       )
		text += 'outputEnergies                  {:<10}                                \n'   .format(self.outputEnergies)
		text += 'outputPressure                  {:<10}                                \n' .format(self.outputPressure)
		
		
		text += '\n\n'
		text += '#############################################################                \n'
		text += '## EXECUTION SCRIPT                                        ##                \n'
		text += '#############################################################               \n\n'
		
		
		if  minimize > 0:
			text += '# Minimization                                                           \n'
			text += 'minimize                        {:<30}  \n'.format(minimize)
			#       'amber                           yes    \n'
		if  dynamics > 0 :
			text += '# Molecular Dynamics                                                     \n'
			#       'amber                           yes    \n'
			text += 'run                             {:<30}  \n'.format(dynamics)

		#print text
		configFile.write(text)


	def set_namd_path (self, namdpath = None, CPUs = None, GPUs =  None):
		""" Function doc """
		self.namd_path  = namdpath #'/home/peagah/NAMD_2.14_Linux-x86_64-multicore/namd2'
		self.nCPUs      = CPUs     #' +p 12 '
		self.nGPUs      = GPUs     #' +devices 0'


	def get_cell (self, inputPDB =  None, computePME = True, verbose = True):
		""" Function doc """

		if inputPDB == None:
			filein = open(self.coordinates, 'r')
			coordType = self.coordinates.split('.')
		else:
			filein = open(inputPDB, 'r')
			coordType = inputPDB.split('.')
		

		#defining coordinate type
		#coordType = filein.split('.')
		#print len(filein)
		if coordType[-1] == 'crd' or coordType[-1] == 'inpcrd' or coordType[-1] == 'rst':
			coords_x = []
			coords_y = []
			coords_z = []
			lines = filein.readlines()
			print(len(lines))
			lines.pop(-1)
			#lines.pop(0)
			for line in lines:
				
				line2 = line.split()
				print(line2)

				if len(line2) == 6:
			
					try:
						#print line[30:38], line[38:46], line[46:54]
						x = float(line2[0])
						y = float(line2[1])
						z = float(line2[2])
						coords_x.append(x)
						coords_y.append(y)
						coords_z.append(z)
						
						x = float(line2[3])
						y = float(line2[4])
						z = float(line2[5])
						coords_x.append(x)
						coords_y.append(y)
						coords_z.append(z)
					except:
						pass
		
		if coordType[-1] == 'pdb':
			coords_x = []
			coords_y = []
			coords_z = []

			for line in filein:
				line2 = line.split()
				#print line2
				try:
					if line[0:4] == 'ATOM' or line[0:4] == 'HETA':
						print(line[30:38], line[38:46], line[46:54])
						x = float(line[30:38])
						y = float(line[38:46])
						z = float(line[46:54])
						
						coords_x.append(x)
						coords_y.append(y)
						coords_z.append(z)
						#print line2, x, y, z
					else:
						pass
				except:
					pass

		coords_x.sort()
		coords_y.sort()
		coords_z.sort()

		#print coords_x[0],coords_x[-1], coords_x[-1]-coords_x[0] 
		#print coords_y[0],coords_y[-1], coords_y[-1]-coords_y[0] 
		#print coords_z[0],coords_z[-1], coords_z[-1]-coords_z[0] 


		self.cellBasisVector1  = [ coords_x[-1]-coords_x[0],                      0. ,   0.                       ]
		self.cellBasisVector2  = [ 0.                      ,coords_y[-1]-coords_y[0] ,   0.                       ]
		self.cellBasisVector3  = [ 0.                      ,                 0.      ,   coords_z[-1]-coords_z[0] ]
		self.cellOrigin        = [ (coords_x[-1])/2        , (coords_y[-1] )/2       ,  (coords_z[-1] )/2         ]
		
		if computePME == True:
			self.PMEGridSizeX      = int((coords_x[-1]-coords_x[0])+1)
			self.PMEGridSizeY      = int((coords_y[-1]-coords_y[0])+1)
			self.PMEGridSizeZ      = int((coords_z[-1]-coords_z[0])+1)
	
		if verbose == True:
			print('cellBasisVector1  {:10}  {:10}  {:10} '.format(self.cellBasisVector1[0], self.cellBasisVector1[1], self.cellBasisVector1[2]))
			print('cellBasisVector2  {:10}  {:10}  {:10} '.format(self.cellBasisVector2[0], self.cellBasisVector2[1], self.cellBasisVector2[2]))			
			print('cellBasisVector3  {:10}  {:10}  {:10} '.format(self.cellBasisVector3[0], self.cellBasisVector3[1], self.cellBasisVector3[2]))			
			print('\ncellOrigin        {:10}  {:10}  {:10} \n'.format(self.cellOrigin      [0], self.cellOrigin      [1], self.cellOrigin      [2]))			
			if computePME == True:			
				print('PMEGridSizeX 		', self.PMEGridSizeX) 			
				print('PMEGridSizeY 		', self.PMEGridSizeY) 			
				print('PMEGridSizeZ 		', self.PMEGridSizeZ) 			
			

	def run_NaMD (self, inputFile):
		""" Function doc """
		
		#self.namd_path  = '/home/peagah/NAMD_2.14_Linux-x86_64-multicore/namd2'
		#self.nCPUs      = ' +p 12 '
		#self.nGPUs      = ' +devices 0 '
		
		cmd = namd_path + " +p 12 " + str(self.CPUs) +" +devices " + str(self.nGPUs) + " " + inputFile + ' > ' + inputFile +'.log'
		os.system(cmd)
