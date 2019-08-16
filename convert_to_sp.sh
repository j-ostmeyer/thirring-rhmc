#!/bin/bash
MODULE=$1
MODULE_SP=${MODULE%.F90}\_sp.F90

cp $MODULE $MODULE_SP

SUBROUTINES=$(cat $MODULE_SP | grep -o -E 'subroutine\s+\w+\s?\(' | sed -E 's/subroutine\s+(\w+)\s?\(/\1/')

KNOWN_TARGET=NO

if  grep "dirac.F90" <(echo $MODULE) &> /dev/null 
then
# Neglected: converting numerical constants.
  KNOWN_TARGET=YES
  sed -i -E 's/module\s+(\w+)\s?$/module \1_sp/;
    s/subroutine\s+(\w+)\s?\(/subroutine \1_sp(/;
    s/end\s+subroutine\s+(\w+)/end subroutine \1_sp/;
    s/complex\s*\(\s*dp\s*\)/complex(sp)/g;
    s/real\s*\(\s*dp\s*\)/real(sp)/g;' $MODULE_SP
elif   grep "comms4.F90" <(echo $MODULE) &> /dev/null  ||  grep "comms5.F90" <(echo $MODULE) &> /dev/null 
then
  KNOWN_TARGET=YES
  sed -i -E 's/module\s+(\w+)\s?$/module \1_sp/;
    s/complex\s*\(\s*dp\s*\)/complex(sp)/g;
    s/real\s*\(\s*dp\s*\)/real(sp)/g;
    s/MPI_Double_Complex/MPI_Complex/g;
    s/(halo_[456]_[xyzt](up|dn)_(send|recv))\(/\1_sp(/g' $MODULE_SP
fi

for SUBROUTINE in $SUBROUTINES
do 
   sed -i 's/\b'$SUBROUTINE'\b/'$SUBROUTINE'_sp/g' $MODULE_SP
done 

if [ "$KNOWN_TARGET" == NO ] 
then
  echo "Error, target not recognized"
  exit 1
fi
# NEGLECTING UNCOMMS_[456]
