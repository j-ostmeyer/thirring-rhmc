#!/bin/bash
MODULE=$1
MODULE_SP=${MODULE%.F90}\_sp.F90

cp $MODULE $MODULE_SP

if [ "$MODULE" == "dirac.F90" ]
then
# Neglected: converting numerical constants.
sed -i -E 's/module\s+(\w+)\s?$/module \1_sp/;
    s/subroutine\s+(\w+)\s?\(/subroutine \1_sp(/;
    s/end\s+subroutine\s+(\w+)/end subroutine \1_sp/;
    s/complex\s*\(\s*dp\s*\)/complex(sp)/g;
    s/real\s*\(\s*dp\s*\)/real(sp)/g;' $MODULE_SP
elif [ "$MODULE" == "comms5.F90" ] 
then
sed -i -E 's/module\s+(\w+)\s?$/module \1_sp/;
    s/subroutine\s+(\w+)\s?\(/subroutine \1_sp(/;
    s/end\s+subroutine\s+(\w+)/end subroutine \1_sp/;
    s/complex\s*\(\s*dp\s*\)/complex(sp)/g;
    s/real\s*\(\s*dp\s*\)/real(sp)/g;
    s/MPI_Double_Complex/MPI_Complex/g;
    s/(halo_5_[xyzt](up|dn)_(send|recv))\(/\1_sp(/g' $MODULE_SP
fi

