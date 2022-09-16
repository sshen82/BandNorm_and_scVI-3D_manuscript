#!/bin/bash

FILES=`ls /p/keles/schic/volumeA/Script/insulation/3DVI/*.txt`
export PERL5LIB=${PERL5LIB}:/u/s/s/sshen82/Rfile/Perl/lib/perl5
export PATH=/u/s/s/sshen82/bedtools2/bin:$PATH

for f in $FILES
do
    perl /u/s/s/sshen82/Rfile/cworld/scripts/perl/matrix2insulation.pl -i $f -o $f -v --is 10000000 --ids 1000000
done
