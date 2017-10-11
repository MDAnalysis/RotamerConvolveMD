#!/bin/bash
# Run the PepTSO example and compare to reference result.
# NOTE: This should be a Python unit test!

rm -rf dcd dat
mkdir dcd dat || { echo "Failed to create work dirs"; exit 1; }
if [ "`which convolve-mtss-rotamers.py`" = "" ]; then
    echo "Could not find convolve-mtss-rotamers.py."
    echo "Install the package first with 'python setip.py install --user' from the top dir"
    exit 2
fi

convolve-mtss-rotamers.py \
    --resid 47 330  \
    --histogramBins 0 80 1  \
    --clashDistance 2.2  \
    --output "dat/peptso-xrd" \
    --plotname "dat/peptso-xrd.pdf" \
    --outputRawDistances "dat/peptso-xrd" \
    --dcdfilename "dcd/peptso-xrd" \
    --dcdfilenameNoClashes "dcd/peptso-xrd" \
    --useNOelectron \
    --libname "MTSSL 298K 2015" \
    peptso.gro 


diff reference_NO_2015/peptso-xrd-47-330.dat dat/peptso-xrd-47-330.dat
test $? -eq 0 && echo "Test PASSED" || echo "Test FAILED."

