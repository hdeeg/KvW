This directory contains a demonstrator
of KvW.pro to determine min-times of two eclipses of CM Dra
 in TESS light curves, with the modified KvW method by Deeg 2020.
 The two light curves are:
 - Primary eclipse at epoch 7024, generating also Fig 1 and first entry
    in Table 1  
 - Incomplete primary at epoch 7023, generating also Figs.3 and 4
 Printed output gives also the timing error from the original equation by 
 Kvw 1956  (with a NaN at epoch 7024)

To execute:
IDL> .rnew kvwdemo1

Dependencies:
rdtab.pro    (generic ASCII table reader)
