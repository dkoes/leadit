leadit
=========

This is *alpha* software. Many features are missing or broken.

The intent for leadit is that it will enable lead optimization of existing
hits by deconstructing a ligand's chemistry and, using the same reaction pathway,
reconstruct alternative leads that match specified shape and pharmacophore criteria.

leadit constructs databases of molecular fragments.  Each database is
constructed around a single reaction framework.  It is populated using
conformers that are compatible with that framework.  A minimal scaffold is
identified from the reaction.  Scaffold conformers are clustered with respect
to RMSD and connection point deviations.  Reaction components are clustered
with respect to unminimized RMSD after alignment to the appropriate scaffold.
 Only representative fragments are retained.

The fragment conformations will be indexed using shapedb with some color
 information.  I may also add in a kd-tree lookup if the color shapedb color
 scheme doesn't work as well as I hope.

 It will also be possible to print statistics of the current database
 and to deconstruct single compounds.
 
 Search functionality for single and multiple databases will be implemented.
 
  I will need to add robust support for I/O parallelism.
 
  Eventually a server mode will be added.

# Installation

An (incomplete) list of dependencies:

`apt install libglpk-dev liblemon-dev libann-dev coinor-*`

RDKit is required.  However, I think a certain amount of effort needs to be
spent verifying and fixing RDKit's reaction handling still.  Also, modern RDKit
doesn't work - have to use my current fork.

OpenBabel is required, but it needs to be the version in the development branch.


The LEMON Graph Library is required (needs to be built from source): http://lemon.cs.elte.hu/trac/lemon/wiki/Downloads

To build use cmake version 3.10 or later:
```
mkdir build
cd build
cmake ..
make -j12
```

# Usage
```
USAGE: leadit [options]

OPTIONS:
  Operation to perform:
    -CreateDatabase        - Create a new reaction database
    -AddMolecules          - Add conformers to database (regenerates indices)
    -SearchDatabase        - Search database for leadit query
    -DatabaseInfo          - Print database information
    -LigandInfo            - Print decomposition of passed ligand(s)
    -Server                - Start leadit server
  -connect-cutoff=<number> - Maximum distance allowed between connection points
  -dbdir=<string>          - database directory(s)
  -exc-mol=<string>        - Receptor excluded shape constraint.
  -exc-shrink=<number>     - Amount to reduce excluded shape.
  -force                   - Overwrite any existing database
  -help                    - Display available options (--help-hidden for more)
  -in=<string>             - input file(s)
  -inc-mol=<string>        - Ligand minimum volume shape constraint.
  -inc-shrink=<number>     - Amount to reduce ligand minimum shape.
  -out=<string>            - output file
  -pharma=<string>         - Pharmacophore search constraints
  -reactant-rmsd=<number>  - Maximum RMSD for reactants to be merged
  -ref=<string>            - Reference scaffold for search
  -rpos=<int>              - Reactant position to replace for search
  -rxn=<string>            - reaction SMARTS file
  -scaffold-rmsd=<number>  - Maximum RMSD for scaffolds to be merged
  -verbose                 - verbose output
```

## Database Creation

`leadit -CreateDatabase  -dbdir db -rxn r.smart -in mols.sdf`

Note that this is quite fragile.  If rdkit fails to decompose a molecule it dies.
It will create a searchable database in the `db` directory.

## Database Info

```
leadit -DatabaseInfo -dbdir db
Reaction: [#6:3]-[#7:1]-[#6:2](=[O:15])-,:[C:5](-,:[#6:4])(-,:[#1,#6:13])-,:[#7:6]-[#6:7](-[#6:8])-[#6:9](=[O:10])-[#7:16](-,:[C:12])-,:[#1,#6;A:11]>>[#6:8]-[#6:7](-[#7:6])-[#6:9](-[#8:14])=[O:10].[#6:4]-[#6:5](-[#1,#6:13])=[O:15].[#6:3]-,:[N&+:1]#[C&-:2].[#1:17]-,:[#7:16](-,:[C:12])-,:[#1,#6;A:11]
Product: [NH:1]([C:2](=[O:15])[C:5]([H:13])([CH3:4])[NH:6][CH:7]([CH3:8])[C:9](=[O:10])[N:16]([H:11])[CH3:12])[CH3:3]
Core: [O:15]=[CH:2][CH2:5][NH:6][CH2:7][CH2:9][NH2:16]

2 total conformers available
```

```
leadit -LigandInfo -in mols.sdf  -dbdir db
[O:15]=[C:2][C:5][N:6][C:7][C:9][N:16]	CCCCc1cccc[c:8]1.[O:10] C[C:4]C CC1c2cc(O)c(O)cc2C[C:3][N:1]1 CN1CC2CC1CN2C(=O)[C:12]c1ccccc1 
[O:15]=[C:2][C:5][N:6][C:7][C:9][N:16]	C[C:8]C.[O:10] CC[C:13]Cc1cccc[c:4]1 CN1CC2CC1CN2C(=O)[C:3](c1ccccc1)[N:1] C[C:11]c1cc(O)c(O)cc1C[C:12] 
```

## Database Search

_This may or may not work.  Pharamacophore constraints definitely aren't implemented_.

```
leadit -SearchDatabase -in one.sdf  -dbdir db -ref one.sdf -rpos 1 -inc-mol one.sdf -inc-shrink 2 -exc-mol rec.pdb -exc-shrink 2
```

This will search the database `db` for modification of the reference moleclue (`-ref one.sdf`) at position 1 in the reaction.  
It will identify all the the modifications at that position that are sterically between the inclusive molecule (with its volume
shrunk by 2) and the recetor (which its volume shrunk by 2).  It's not clear this is currently working (in particular the steric constraints).

If you specify an output file, it will output the hits, but only the fragment part and not in a particularly useful location. :-(



# License

GNU GENERAL PUBLIC LICENSE Version 2.
