#include "arguments.h"
#include <iostream>

int arguments::help_print() {
  if (help == false)
    return 1;
  else {
    cout << R"(
================================================================================
                       TITAN ENGINE — USER MANUAL
================================================================================

SYNOPSIS
    elec_vdw.exe [OPTIONS]
    elec_vdw.exe -a arg.in

DESCRIPTION
    Titan is a high-performance molecular analysis engine for computing
    electrostatic, van der Waals interactions, GBSA solvation energies,
    hydrogen bonds, dynamic cross-correlation, PCA, TICA, MSM, LLPS analysis
    and more from PDB/XTC/TRR trajectory files.

================================================================================
1. INPUT / OUTPUT
================================================================================

  -f <file>             Single PDB/GRO input file (topology reference).
                        Default: complex.pdb

  -i <file>             List file containing multiple trajectory filenames.
                        Default: input.txt
                        Each line = one trajectory file (PDB or XTC/TRR).

  -traj <file>          XTC or TRR binary trajectory file (Gromacs format).
                        Auto-detects format from extension (.xtc or .trr).
                        Requires -f for topology reference.
                        Coordinates: nm -> Angstrom (x10) auto-converted.
                        arg.in syntax: -traj=md.xtc

  -1pdb                 Use single-PDB mode (all frames in one file).

  -a <file>             Read arguments from file (default: arg.in).
                        All options below can be placed in this file.
                        Note: in arg.in, replace ' ' by '='.
                        Example: -traj=md.xtc, -mol2=7589

  -opdb [name.pdb]      Output PDB file. Default: out_complex.pdb
  -opdbsep [name.pdb]   Output PDB with separated results per trajectory.

  -n <file>             Index file. Default: index.ndx

  -topo <file>          Topology file for molecule 1.
                        Default: lig.itp. Supports .itp, .top, .parm7

  -topo2 <file>         Topology file for molecule 2. Default: no

  -ff <name>            Force field name. Default: AMBER99SB
                        Supported: AMBER99SB, AMBER99SB-ILDN, GROMOS43A1,
                        GROMOS53A6, GROMOS54A7, CHARMM27, CHARMM36,
                        OPLSAA, MARTINI

  -pattern-input <pat> [b] [e] [step]
                        Pattern input mode. <pat> must contain *.
                        Tool replaces * with numbers from b to e.
                        Default: b=0, step=1.

================================================================================
2. GROUPS & MOLECULES
================================================================================

  -g <G1>&<G2>          Interaction groups (e.g. -g 1&2 or -g 1-2).
  -sg <SG1>-<SG2>       Sub-groups within main groups.

  -mol-1-type <type>    Type of molecule 1.
                        Options: protein, complex, lig, rna, dna
                        Default: protein

  -mol-2-type <type>    Type of molecule 2.
                        Options: protein, protein-user, complex, lig, rna, dna
                        Default: lig

  -mol2 <N>             First atom number of molecule 2.

  -n_wat <N>            Number of water molecules. Default: -1 (off)

================================================================================
3. GBSA / SOLVATION
================================================================================

  -gbsa                 Enable GBSA solvation energy calculation.
  -alpb                 Enable Analytical Linearized PB (implies -gbsa).
  -igb <N>              GB model. Default: 5
  -ext_diel <val>       External dielectric constant. Default: 78.5
  -int_diel <val>       Internal dielectric constant. Default: 1.0
  -saltcon <val>        Salt concentration (M). Default: 0.15
  -epsilon_in <val>     Epsilon_in parameter. Default: 1.0

================================================================================
4. PERIODIC BOUNDARY CONDITIONS (PBC)
================================================================================

  -pbc                  Enable basic PBC mode.
  -pbci <chain_ref>     PBC with internal reference chain.
  -pbci2 <chain_ref>    PBC mode 2 (more aggressive).
  -pbci3 <chain_ref>    PBC mode 3 (most aggressive).
  -pbca <chain_ref>     Aggressive PBC mode.
  -pbcn                 PBC normal mode.
  -pbcclust <axes>      PBC clustering mode. <axes> = xyz combination.

================================================================================
5. ENERGY & FORCE MAPS
================================================================================

  -emap [type]          Enable energy map calculation.
                        Types: residue-atom (default), residue-residue,
                               atom-atom, atom-residue
  -emap-level <level>   Energy map level: per-res (default) or all.

  -fmap [type]          Enable force map calculation. Same types as -emap.
  -fmap-level <level>   Force map level: per-res (default) or all.

  -energy-res           Print per-residue energy decomposition.
  -force-res            Print per-residue force decomposition.
  -new-emap             Use new energy map format.

  -com_energy <axis>    COM-distance energy projection.
                        Options: x, y, z, all (distance), or default.

  -denfo <val>          Distance-energy cutoff (Angstrom). Default: infinity
  -dsc <val>            Standard distance cutoff. Default: 6.5 Angstrom

================================================================================
6. ANALYSIS MODULES
================================================================================

  -network              Enable DCCM network analysis.
  -pca                  Enable PCA (implies -network).
  --fes-grid=<N>        FES grid resolution. Default: 60

  -allostery            Enable Information-Theory allosteric analysis
                        (Generalized Correlation + Transfer Entropy).
                        Implies -network.
  --allostery-tau=<N>   TE time lag in frames. Default: 5
  --allostery-bins=<N>  Histogram bins (-1 = adaptive). Default: -1

  -tica                 Enable TICA (implies -network).
  --tica-lag=<N>        TICA lag time (frames). Default: 10
  --tica-components=<N> TICA components. Default: 5

  -msm                  Enable MSM (implies -tica, -network).
  -msm-clusters <N>     K-means clusters. Default: 100
  -msm-lag <N>          MSM lag time (frames). Default: 10
  -msm-nts <N>          Number of implied timescales. Default: 3
  -msm-kmeans-iter <N>  K-means iterations. Default: 300
  -msm-seed <N>         Random seed. Default: 42
  -msm-macrostates <N>  Macrostates (-1 = auto). Default: -1
  -msm-no-boundary-split  Disable boundary splitting.
  -msm-no-pcca          Disable PCCA analysis.
  -msm-no-its           Disable implied timescales.
  -msm-its-lags <list>  ITS lag list (comma-separated, e.g. 1,2,5,10).

  -llps                 Enable LLPS phase separation analysis.
  -llps-rcut <val>      Contact cutoff (nm). Default: 0.7
  -llps-perc <val>      Percolation threshold. Default: 0.5
  -llps-maxlag <N>      Max lag frames for C(t). Default: 500
  -llps-dimer-size <N>  Chains per macro-node. Default: 1
  -llps-contacts        Enable residue contact map output.
  -llps-sticker <beads> Comma-separated bead names.
                        Default: BB,SC1,SC2,BB1,BB2,BB3,SC3,SC4

  -cas                  CAS analysis.
  -tda                  TDA analysis.
  -ml                   ML affinity prediction.

================================================================================
7. HYDROGEN BONDS
================================================================================

  -hb                   Enable H-bond calculation.
  -hb side-chain        H-bond for side-chain only.
  -ahb <val>            H-bond angle cutoff (degrees). Default: 135.0
  -dhb <val>            H-bond distance cutoff (Angstrom). Default: 3.5

================================================================================
8. ENTROPY
================================================================================

  -entropy              Enable interaction entropy (Duan-Zhang).
  -entropy-quasi <com> <lig> <rec>
                        Quasi-harmonic entropy with group indices.

================================================================================
9. TRAJECTORY / TIME CONTROL
================================================================================

  -b <time>             Begin time (ps). Default: 0
  -e <time>             End time (ps). Default: very large
  -dt <N>               Frame step interval. Default: 0 (use all)
  -nt <N>               Number of threads. Default: auto

================================================================================
10. FORCE FIELD / PARAMETERS
================================================================================

  -path-a <dir>         Path to amino acid parameters. Default: lib_ff/
  -path-n <dir>         Path to nonbonded parameters. Default: lib_ff/
  -eu <unit>            Energy unit. Default: kJ/mol
  -fu <unit>            Force unit. Default: pN

================================================================================
11. MISCELLANEOUS
================================================================================

  -nocal                No calculation mode (parse/PBC/HB only for 1 PDB).
  -cscd                 Calculate side-chain contact distances.
  -pqr                  Output PQR format.
  -cal-force            Calculate atomic forces.
  -cal-volume           Calculate molecular volume.
  -max-force            Report maximum force (auto-enables -cal-force).
  -cache / -cache -no   Enable / disable cache mode.
  -use-CI               Use 95% confidence intervals instead of std dev.
  -elec-screening       Enable electrostatic screening.
  -interchain-energy    Calculate inter-chain energies.
  -dipole-moment        Calculate dipole moments.
  -split                Split results per trajectory.
  -all-subfolders       Process all sub-folders (implies -split).
  -cat-result           Concatenate results from multi-input files.
  -all-respair          Output all residue pairs.
  -resdis2com           Residue distance to COM.
  -pull-direct          Find pulling direction.
  -charged_atom_dis [N] Charged atom distance distribution.
  -water-type <type>    Water analysis type: surface or sphere.
  -dis-water <val>      Water distance cutoff. Default: 3.5

  -h, --help            Print this help message.
  -v, --version         Print version number.

================================================================================
12. EXAMPLES
================================================================================

  # Basic energy decomposition (PDB):
  elec_vdw.exe -f complex.pdb -g 1&2 -ff AMBER99SB

  # GBSA solvation energy:
  elec_vdw.exe -f complex.pdb -g 1&2 -gbsa -igb 5

  # XTC trajectory with PDB topology:
  elec_vdw.exe -f reference.pdb -traj md.xtc -g 1&2

  # Multi-file XTC (via input.txt):
  elec_vdw.exe -f reference.pdb -traj md_1.xtc -i input.txt -g 1&2

  # LLPS phase separation (Martini):
  elec_vdw.exe -f system.pdb -traj md.xtc -llps -ff MARTINI

  # Full analysis pipeline:
  elec_vdw.exe -f complex.pdb -g 1&2 -gbsa -network -pca -tica -msm

  # Using arg.in file:
  elec_vdw.exe -a arg.in

  # Example arg.in:
  #   -f=complex.pdb
  #   -g=1&2
  #   -ff=AMBER99SB
  #   -gbsa
  #   -network
  #   -traj=md.xtc
  #   -llps
  #   -nt=8

================================================================================
)";
    return 0;
  }
}
