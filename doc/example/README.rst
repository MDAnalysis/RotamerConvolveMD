Running the PepTso example
==========================

DEER
----

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
        --libname "MTSSL 298K 2015" \
        peptso.gro 

Compare the output to the reference::

   diff ../../tests/mtssl_2015/reference_NO_2015/peptso-xrd-47-330.dat  dat/peptso-xrd-47-330.dat

(If there is no difference then no output will be shown --- this is
the expected result.)   

You can also look at the PDF of the distance histogram and compare to
``reference/peptso-xrd-47-330.pdf``.

The above test can be run by executing the script
``tests/mtssl_2015/test_mtssl_2015.sh`` from the top directory.


PRE
---

Run the PRE example from this directory::

   mkdir dcd dat
   convolve-mtss-rotamers_pre.py \
        --resid 47  \
        --clashDistance 2.2  \
        --plotname "dat/peptso-xrd.pdf" \
        --outputRawDistances "dat/peptso-xrd" \
        --dcdfilenameAll "dcd/peptso-xrd" \
        --dcdfilenameNoClashes "dcd/peptso-xrd" \
        --libname "MTSSL 298K 2015" \
        peptso.gro 

Compare the output to the reference::

   diff ../../tests/mtssl_2015/reference_pre_NO_2015/peptso-xrd-47-rawDistances.dat dat/peptso-xrd-47-rawDistances.dat

(If there is no difference then no output will be shown --- this is
the expected result.)   
   
You can also look at the PDF of the distance histogram and compare to
``reference_pre/peptso-xrd-47.pdf``.

The above test can be run by executing the script
``tests/mtssl_2015/test_mtssl_pre_2015.sh`` from the top directory.

  
