# The 27 proteins of the paper

The 27 proteins of Figure 5 and Supplementary Table S7 of the paper, the benchmark set of the FAGBI study, from
the deposited structures to the solvation energies.

| File | Contents |
|---|---|
| `inputs/<id>_source.pdb` | The structure as deposited in the PDB (https://files.rcsb.org/download/<ID>.pdb) |
| `inputs/<id>_model1_all_protein.pdb` | Its first model, `ATOM` records only |
| `inputs/<id>_full_charmm_ph7.pqr` | The PDB2PQR 3.7.1 output used by the paper: hydrogens added, PROPKA protonation at pH 7.0, CHARMM charges and radii |
| `inputs/<id>_pdb2pqr.log` | The PDB2PQR log |
| `inputs/fagbi27_fixture_manifest.csv` | Source URL, preparation policy, PDB2PQR options and SHA-256 of every file |
| `prepare_pqr.py` | The preparation (section 3 of `MANUAL_BIE_PB.md`); `--verify` rebuilds the 27 PQR files from the deposited PDB files and checks them byte by byte |
| `protocol_arg.in.template` | The argument file of every protein calculation of the paper; `{PQR}` stands for the input |
| `run_all.py` | Runs the protocol on all 27 proteins and compares with the paper |
| `reference_protocol.csv` | The paper's results: node patches, GMRES iterations, ΔG_solv (kcal/mol), true residual and sampled operator discrepancy e_op,S |

1CBN has alternate conformers that PDB2PQR rejects; it was prepared with `prepare_pqr.py --policy clean` (input
hydrogens removed, one conformer per residue ranked by occupancy). The other 26 proteins keep their `ATOM`
records as deposited.

```bash
pip install pdb2pqr==3.7.1             # only for prepare_pqr.py
python prepare_pqr.py --verify         # optional: rebuild the 27 PQR files
python run_all.py --titan ../../../titan.exe
```

`run_all.py` writes one directory per protein under `run/` and the table `run/results.csv`. With this release
the 27 energies, iteration counts and node-patch counts are those of the paper. 451C does not pass the operator
gate (e_op,S ≤ 8 × 10^-3) at h = 1.25 Å; the paper reports it with this failure and repeats it at h = 0.75 and
0.5 Å, where the gate passes (Section 3.4).
