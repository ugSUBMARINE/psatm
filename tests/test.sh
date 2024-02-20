#!/bin/bash

nomiss=$(../target/debug/psatm ../example/gb1.pdb ./gb1_test_out.pdb)
onemiss=$(../target/debug/psatm ../example/gb1_missing.pdb ./gb1_missing_test_out.pdb )

if [[ ! $(echo $nomiss | wc -l) -eq 1 ]]; then
    echo "FAILED: Warning when could not generate all pseudoatoms even no atoms are missing"
    echo $nomiss
    exit 1
else
    echo "PASS: Warning pseudoatom with no missing atoms"
fi

if echo $onemiss | grep -q 'Not able to calculate pseudoatom for \[ GLU A  42 \]'; then
    echo "FAILED: Warning for wrong output for one missing side chain (GLU A 42) that should yield in one impossible calculation"
    echo $onemiss
    exit 1
else
    echo "PASS: Warning pseudoatom calculation with one missing side chain"
fi

if ! cmp --silent ../example/gb1_out.pdb ./gb1_test_out.pdb;then
    echo "FAILED: no missing sidechain generates wrong output file"
    exit 1
else
    echo "PASS: No missing side chain file generated"
fi

if ! cmp --silent ../example/gb1_missing_out.pdb ./gb1_missing_test_out.pdb;then
    echo "FAILED: missing sidechain generates wrong output file"
    exit 1
else
    echo "PASS: missing side chain file generated"
fi
