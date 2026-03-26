#include<iostream>
#include"arguments.h"
int arguments::help_print()
{
	if (help == false)
		return 1;
	else
	{
		cout << "This tool uses .PDB file. You can use 1 .pdb file or 2 .pdb files.\nIf you use only 1 .pdb file, you can specify it by -f option\n";
		cout << "If you use 2 .pdb files, your .PDB input files must have format: lg_xxx.pdb and rec_xxx.pdb," << '\n';
		cout << " in which lg is ligand (or protein 1st) and rec is receptor (or protein 2nd)" << '\n';
		cout << "xxx is index of snaps." << '\n';
		cout << "If you want to use input.txt, your input.txt must have format: each index is written line by line." << '\n';
		cout << " If you do not use input.txt, your snap index must begin by 0." << '\n';
		cout << "+ This application also supports other features, please consult the usage\n";
		cout << "Topology file for ligand (molecule 1) is group 2, topology file for receptor (molecule 1) is group 1.\n";
		cout << "*************************************************************************************************************" << '\n';
		cout << "Note: Your system should have multi-cores processor for parallel calculation, especially with large molecules\n";
		cout << "*************************************************************************************************************" << '\n';
		cout << "USAGE: -v to print software build number\n"
			<< "       -ff <forcefield> (supported forcefields: AMBER99SB, AMBER99SB-ILDN, GROMOS43A1, GROMOS53A6, GROMOS54A7, CHARMM27, CHARMM36 and OPLSAA)\n"
			<< "       -cscd to monitor distance between com of all residues of receptor and com of all residues/atoms of ligand\n"
			<< "       -all-subfolders to use all .pdb files in folder and its subfolders\n"
			<< "       -split to print result to folder(s) with folder(s)'s name is .pdb file name\n"
			<< "       -pqr to create .PQR files\n"
			<< "       -act to calculate all side-chain contact\n"
			<< "       -pbc to activate built-in pbc solver\n"
			<< "       -pbci to activate intensive built-in pbc solver\n"
			<< "       -pbci <chain id> to activate new intensive built-in pbc solver, chain id begin from 0 (0 means chain A)\n"
			<< "       -pbca to activate aggressive pbc solver \n"
			<< "       -pbc-envi to specify chains which are enviroment (do not try to solve pbc for these chains just residues), useful for water, lipid molecules\n"
			<< "       -1pdb to use only 1 pdb file for input instead of lg_*.pdb and rec_*.pdb\n"
			<< "       -charged_atom_dis <binnum> to calculate histogram of distances between charged atoms and center of mass of molecule\n"
			<< "       -useForceOptionforEnergy to use force options for energy options\n"
			<< "       -cal-force to calculate interaction force\n"
			<< "       -max-force to calculate maximum force, this option will automatically activate -cal-force\n"
			<< "       -emap to create energy map, or -emap <residue-atom>/<residue-residue>/<atom-residue> or <atom-atom> to specify type\n"
			<< "       -fmap to create force map, or -fmap <residue-atom>/<residue-residue>/<atom-residue> or <atom-atom> to specify type\n"
			<< "       -emap-level <all/per-res>to specify level of emap calculation\n"
			<< "       -fmap-level <all/per-res>to specify level of fmap calculation\n"
			<< "       -new-emap to print average energy interaction for per residue\n"
			<< "       -hb to find hydrogen bonds (HBs), -hb side-chain to find HBs of side chains only\n"
			<< "       -entropy to calculate entropy by duan-zhang method\n"
			<< "       -nocal to generate .pqr, solve pbc problem or find hydrogen bonds for only 1 pdb file\n"
			<< "       -energy-res to print interaction energy of residues, note: it will enable -emap\n"
			<< "       -force-res to print interaction force of residues, note: it will enable -fmap\n"
			<< "       -com_energy <axis> <axis=x/y/z/all> to print interaction energy between coms of 2 mols along <axis> (all for distance)\n"
			<< "       -side_chain_contact full-res to calculate side chain contact for whole residue\n"
			<< "       -opdb <filename> to write pdb file of system\n"
			<< "       -use-CI to calculate confidence interval 95% insteads of standard deviation\n"
			<< "       -a <file name> is argument file, default is arg.in\n"
			<< "       -i <file name> is pdb list file, default is input.txt\n"
			<< "       -b <time (ps)> first frame to read from input file\n"
			<< "       -e <time (ps)> last frame to read from input file\n"
			<< "       -dt <time (ps)> is time interval between snapshots (frame)\n"
			<< "       -eu <energy unit> is energy unit, this tool supports kJ/mol and kcal/mol, default is kJ/mol\n"
			<< "       -fu <force unit> is force unit, this tool supports kJ/mol/nm and pN, default is pN\n"
			<< "       -g <group 1>&<group 2> or <group 1>-<group_2> to choose groups to calculate\n"
			<< "       -sg <mol1-subgroup1> <mol1-subgroupN> - <mol2-subgroup1> <mol2-subgroupN> to calculate energy, force for subgroups\n"
			<< "       -nt <number of threads> to choose number of threads\n"
			<< "       -mol2 <atom number> is the first atom number of molecule 2\n"
			<< "       -topo <file name> is ligand topo file, default is lig.itp. This program supports .itp, .top and .parm7\n"
			<< "       -topo2 <file name> is molecule 2 topo file. This program supports .itp, .top and .parm7\n"
			<< "       -f <file name> is your pdb input file, default is complex.pdb\n"
			<< "       -n <file name> is your index file, default is index.ndx\n"
			<< "       -dsc <distance (Angstrom)> is distance cutoff for side-chain contact calculation, default is 6.5 Angstrom\n"
			<< "       -denfo <distance (Angstrom)> is distance cutoff for energy/force calculation, default is infinity\n"
			<< "       -path-a <path> is the path to aminoacids parameter files\n"
			<< "       -path-n <path> is the path to nonbonded parameter files\n"
			<< "       -ahb <angle (degree)> minimum angle of donor-hydrogen-aceptor, default is 135 degree\n"
			<< "       -dhb <distance (Angstrom)> maximum distance between donor-aceptor, default is 3.5 Angstrom\n"
			<< "       -entropy-quasi <group 1> <group 2> <group 3> to calculate entropy by quasi harmonic method (experimental-for testing only)\n"
			<< "       -pattern-input <pattern> <begin> <end> <step> to calculate multi input file with same pattern\n Your pattern must contain *, this tool will replace * by numbers from begin to end\n"
			<< "If begin, end and step are not provided, default value of begin is 0 and step is 1\n"
			<< "Note: if you use argument file, replace ' ' by =, for ex.: -a=arg.in -pattern-input=com*.pdb 0 10 2, -mol2=7589, -entropy-quasi=2 3 4\n";
		return 0;
	}
}