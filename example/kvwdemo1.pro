;pro KVWdemo1
; Demonstrator of KvW.pro to determine min-times of two eclipses of CM Dra
; in TESS light curves, with revised KvW method by Deeg 2020.
; The two light curves are:
; - Primary eclipse at epoch 7024, generating also Fig 1 and first entry
;    in Table 1  
; - Incomplete primary at epoch 7023, generating also Figs.3 and 4
; Printed output gives also the timing error from the orignal equation by 
; Kvw 1956  (with a NaN at epoch 7024)
;HJD 17nov2020, 

infile =strarr(2)
infile[0]='CMDra7024.lc'
infile[1]='CMDra7023.lc'
;both light curves were
;extracted from PDCSAP_FLUX of MAST file
;tess2019253231442-s0016-0000000199574208-0152-s_lc.fits
;and processed as described in the paper

for lc=0,1 do begin
   rdtab,infile[lc],time,flux,comc='#'   ;invokes ascii table reader
;   p1=plot(time,flux,'-',xtitle='BJD',ytitle='flux',title=infile[lc])

   rms=0.00138   ;rms of off-eclipse data, determined from revision of all
;  Sector 16  CM Dra eclipses, see paper
   mintim=kvw(time,flux,minerr=minerr,kvwminerr=kvwminerr,rms=rms,nfold=5,init_minflux=1,noprint=1,noplot=0)
;init_minflux=1 (using the point of minimum flux  as intial guess) is used 
; above as else the example on the incomplete eclipse would fail.

   print,'infile= ',infile[lc]
   print,format='(a,f15.7,a,f9.7,a,f9.7)','mintime: ',mintim,'+-',minerr,' orig. KvW error: ',kvwminerr
   print,'----------------------------------'
   tmp=""
   if lc lt 1 then read,'hit return for next lightcurve',tmp
endfor  ;lc
end
