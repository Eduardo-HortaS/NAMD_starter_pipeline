## MMPBSA Pipeline Tutorial


# leaprc content 

Loading the forcefield: 

```source PATH/TO/FORCE/FIELD/leaprc.ff99sb```

Followed by loading our pdbs,

Solvating the system: 

```system = loadpdb system.pdb```

Protein in vacuum: 

```protein = loadpdb protein.pdb```

Ligand in vacuum: 

```ligand = loadpdb ligand.pdb```

Complex in vacuum: 

```complex = loadpdb complex.pdb```

Obs: if the ligand was parameterized externally, you'll need to indicate the name of its residues and load their parameters. 

Solvation and neutralization of the system:

```
solvateBox system <box type> <distance between extremities/ends>
addions2 system <ion type> <system charge>
```

Afterwards, you'll have to save your parameters for all charged and treated components:

```
saveamberparm system system.top system.crd
savepdb system system.pdb
saveamberparm protein protein.top protein.crd
saveamberparm ligand ligand.top ligand.crd
saveamberparm complex complex.top complex.crd
quit
```

# build_scripts.py content

Having saved these parameters for all components, it is now time to run the simulation (NAMD) for the system after solvation and neutralization. We're going to use build_scripts.py to generate .namd config files based on the .top, .crd e .pdb files for your system.

```
 crd           = 'system.crd'  
 parm          = 'system.top'
 pdb_fixed1    = 'system.pdb'
```


```
 outputname_00 = 'system1'
 outputname_01 = 'system2'
 outputname_02 = 'system3'
```

Where output 00 corresponds to solvent equilibration, 01 to system equilibration and 02 to production.

# pdb_to_pdb_fixed.py content

System 1 requires water molecule to be free while the rest of the system is fixated for solvent equilibration, making it necessary to perform the following modification of pdb_to_pdb_fixed.py:

```
def export_PDB_fixed (pdbin   = system.pdb,
                      pdbout  = fix.pdb,
                      wat     = True,
                      lig     = False,
                      protein = False,
```

Note that equilibration and production times can be changed.

# NaMD_Object3.py content

Inside NaMD_Object3.py it is possible to change simulation configurations according to user needs. To generate a NAMD config file, use:

```
namd system.namd
```

After the simulation, the resulting system.dcd is a trajectory file that can be used in MMPBSA utilizing its input. 
