# psatm

![TEST](https://github.com/ugSUBMARINE/psatm/actions/workflows/rust.yml/badge.svg)
[![Software License](https://img.shields.io/badge/license-MIT-brightgreen.svg)](/LICENSE.md)
<a title="Code Size" target="_blank" href="https://github.com/ugSUBMARINE/psatm"><img src="https://img.shields.io/github/languages/code-size/ugSUBMARINE/psatm"></a>

Based on a given PDB file, this tool creates a new PDB file that includes pseudoatoms for each residue. These pseudoatoms represent the centroids of each residue's catalytic atoms.

If any atoms are missing, the centroid will be calculated based on the available atoms. If there are no catalytic atoms for a specific residue, a warning will be issued.

![alt text](./example/output.png?raw=true)

## Installation
[Install rust and cargo](https://www.rust-lang.org/tools/install)

Run `cargo build --release` to create a binary of the program.

## Usage
`psatm INFILE OUTFILE`

`psatm example/gb1.pdb example/gb1_out.pdb`

### Examples
Can be found [here](https://github.com/ugSUBMARINE/psatm/tree/master/example).

The respective outputs are the files with the postfix `_out`.

File `gb1_missing.pdb` has the side chain atoms of residue 42 removed to show the output when atoms are missing.
