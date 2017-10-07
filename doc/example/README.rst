Running the PepTso example
--------------------------

Run the DEER example from this directory::

   mkdir dcd dat
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
        --libname "MTSSL 298K 2017" \
        peptso.gro 

Compare the output to the reference::

   diff reference_NO_2017/peptso-xrd-47-330.dat dat/peptso-xrd-47-330.dat

You can also look at the PDF of the distance histogram and compare to
reference/peptso-xrd-47-330.pdf.

The above test can be run by executing the script ``test_mtssl_2017.sh``.


--------------------------


Run the PRE example from this directory::

   mkdir dcd dat
   convolve-mtss-rotamers_pre.py \
        --resid 47  \
        --clashDistance 2.2  \
        --plotname "dat/peptso-xrd.pdf" \
        --outputRawDistances "dat/peptso-xrd" \
        --dcdfilenameAll "dcd/peptso-xrd" \
        --dcdfilenameNoClashes "dcd/peptso-xrd" \
        --libname "MTSSL 298K 2017" \
        peptso.gro 

Compare the output to the reference::

   diff reference_pre_NO_2017/peptso-xrd-47-rawDistances.dat dat/peptso-xrd-47-rawDistances.dat

You can also look at the PDF of the distance histogram and compare to
reference_pre/peptso-xrd-47.pdf.

The above test can be run by executing the script ``test_mtssl_pre_2017.sh``.

  
