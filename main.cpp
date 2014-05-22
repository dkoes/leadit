/*
 * 5/21/2014
 * David Koes (dkoes@pitt.edu)
 *
 * The intent for leadmaker is that it will enable lead optimization of existing
 * hits by deconstructing a ligand's chemistry and, using the same reaction pathway,
 * reconstruct alternative leads that match specified shape and pharmacophore criteria.
 *
 * leadmaker constructs databases of molecular fragments.  Each database is
 * constructed around a single reaction framework.  It is populated using
 * conformers that are compatible with that framework.  A minimal scaffold is
 * identified from the reaction.  Scaffold conformers are clustered with respect
 * to RMSD and connection point deviations.  Reaction components are clustered
 * with respect to unminimized RMSD after alignment to the appropriate scaffold.
 * Only representative fragments are retained.
 *
 * The fragment conformations will be indexed using shapedb with some color
 * information.  I may also add in a kd-tree lookup if the color shapedb color
 * scheme doesn't work as well as I hope.
 *
 * It will also be possible to print statistics of the current database
 * and to deconstruct single compounds.
 *
 * Search functionality for single and multiple databases will be implemented.
 *
 * I will need to add robust support for I/O parallelism.
 *
 * Eventually a server mode will be added.
 */

#include "CommandLine2/CommandLine.h"


enum CommandEnum
{
	create,
	add,
	index,
	search,
	status,
	server
};

cl::opt<CommandEnum> Command(cl::desc("Operation to perform:"), cl::Required,
		cl::values(clEnumVal(create, "Create a new reaction database"),
				clEnumVal(add, "Add conformers to database (invalidates indices)"),
				clEnumVal(index,"Index conformer information"),
				clEnumVal(search,"Search database for leadmaker query"),
				clEnumVal(status,"Print database information"),
				clEnumVal(server,"Start leadmaker server"),
				clEnumValEnd));


int main(int argc, char *argv[])
{
	cl::ParseCommandLineOptions(argc, argv);

}
