#!/bin/bash
# Run the PepTSO example and compare to reference result.
# NOTE: This should be a Python unit test!

rm -rf dcd dat
mkdir dcd dat || { echo "Failed to create work dirs"; exit 1; }
if [ "`which convolve-mtss-rotamers_pre.py`" = "" ]; then
    echo "Could not find convolve-mtss-rotamers_pre.py."
    echo "Install the package first with 'python setip.py install --user' from the top dir"
    exit 2
fi

convolve-mtss-rotamers_pre.py \
    --resid 47  \
    --clashDistance 2.2  \
    --plotname "dat/peptso-xrd-47.pdf" \
    --outputRawDistances "dat/peptso-xrd" \
    --dcdfilenameAll "dcd/peptso-xrd" \
    --dcdfilenameNoClashes "dcd/peptso-xrd" \
    --useNOelectron \
    --libname "MTSSL 298K 2015" \
    peptso.gro 


diff reference_pre_NO_2015/peptso-xrd-47-rawDistances.dat dat/peptso-xrd-47-rawDistances.dat
test $? -eq 0 && echo "Test PASSED" || echo "Test FAILED."
