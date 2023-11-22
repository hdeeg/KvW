;pro KVWdemo1
; Demonstrator of KvW.pro to determine min-times of two eclipses of CM Dra
; in TESS light curves, with revised KvW method by Deeg 2020.
; The two light curves are:
; - Primary eclipse at epoch 7024, generating also Fig 1 and first entry
;    in Table 1  
; - Incomplete primary at epoch 7023, generating also Figs.3 and 4
; Printed output gives also the timing error from the orignal equation by 
; Kvw 1956  (with a NaN at epoch 7024).

;Running this program, with .rnew kvwdemo1,
;the text-output should be:
;infile= CMDra7024.lc
;mintime:   58739.9291169+-0.0000125 orig. KvW error:      -NaN
;----------------------------------
;infile= CMDra7023.lc
;mintime:   58738.6607358+-0.0000191 orig. KvW error: 0.0000662
;----------------------------------


;HJD 17nov2020, first version
;    21nov2023  added expected text-output to comments above.

inpath='../example_data/'
infile =strarr(2)
infile[0]=inpath+'CMDra7024.lc'
infile[1]=inpath+'CMDra7023.lc'
;both light curves were
;extracted from PDCSAP_FLUX of MAST file
;tess2019253231442-s0016-0000000199574208-0152-s_lc.fits
;and processed as described in the paper

for lc=0,1 do begin
   rdtab,infile[lc],time,flux,comc='#'   ;invokes ascii table reader

   rms=0.00138   ;rms of off-eclipse data, determined from revision of all Sector 16  CM Dra eclipses, see paper
   mintim=kvw(time,flux,minerr=minerr,kvwminerr=kvwminerr,rms=rms,nfold=5,init_minflux=1,noprint=1,noplot=0)
;   mintim=kvwcore(time,flux,minerr=minerr,kvwminerr=kvwminerr,rms=rms,nfold=5,init_minflux=1)  ;the same using just the core kvw code

   print,'infile= ',infile[lc]
   print,format='(a,f15.7,a,f9.7,a,f9.7)','mintime: ',mintim,'+-',minerr,' orig. KvW error: ',kvwminerr
   print,'----------------------------------'
   tmp=""
   if lc lt 1 then read,'hit return for next lightcurve',tmp
endfor  ;lc
end
