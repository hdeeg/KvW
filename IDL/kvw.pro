function kvw, time, flux, minerr=minerr, kvwminerr=kvwminerr, init_minflux=init_minflux, rms=rms, nfold=nfold, notimeoff=notimeoff, noplot=noplot, noplnum=noplnum,noprint=noprint, errflag=errflag,debug=debug
  
;Returns eclipse mid (minimum)-time using the Kwee Van-Woerden (1956) method with revised timing error, following Deeg 2021 (Galaxies, vol. 9, issue 1, p. 1)
;Data-points that are equidistant in time are required and the lightcurve should only contain the eclipse (no off-eclipse data).
;For the initial guess of the minimum time, two options are given, controlled by the keyword init_minflux : By default (init_minflux=0), the middle of 
;the lightcurve is assumed. Else (init_minflux=1), the time of the flux-minimum is used. This option is fine or preferential for low-noise lightcurves showing 
;a clear flux-mimimum not extending over more than 2-3 points. 
;If the rms (point-to-point noise) of the off-eclipse flux is known, it should be supplied by keyword rms. Else, rms is estimated from S of the best 
; flux-pairing, which is assumed to be dominated by flux measurement errors. For details on calculation of that rms, see Deeg 2021. 
;By default, the time of minimum is derived from 5 pairings (nfold paramter) that fold at i-1,i-0.5,i,i+0.5,i+1, where i is the index of the point of minimum 
; flux. The original KvW algorithm uses nfold=3.
  
  
;input 
                                ;time   vector with time values in ascending order
                                ;flux   vector with corresponding flux values
;output 
                                ;time of mimimum

;keywords
                                ;minerr  output: error of the minimum time, method by Deeg 2021
                                ;kvwminerr  output: error of the minimum time, original method by KvW 1956
                                ;init_minflux input: By default, the middle of the lightcurve is used as
                                ;        intial guess of the minimum time. If 1, uses as initial min-time estimate the value of lowest flux. 
                                ;rms     input: average measurement error of individual flux values (in units of the flux values). rms should usually be supplied.     
                                ;nfold   input: number of foldings on which to perform pairings of flux values, default=5
                                ;notimeoff input: if set to 1, disables internal offsetting of time values. Is save to use if eclipse min-time is within +-0.1 of time value
                                ;noplot  input: if set to 1, supresses plots
                                ;noplnum input: supresses point-numbers in plot of input lightcurve
                                ;noprint input: if set to 1, supresses all printing
                                ;errflag output flag with
                                ;error-codes. 0: OK, data are equidistant within 1% against medium spacing. 1: Data are not equidistantly spaced
                                ;debug   input debugging flag. If 1, prints some, if 2, prints lots, if 3 even more intermediate values
 
;HJD 19nov2019 initial version
;  24jan2020   introduced internal offset of input times by value of central time-point (floored to first digit) to avoid numerical
;               problems with polyfit on large time-values (even values close to 1 are problematic!). Also new kw notimeoff 
;  31jul2020   added errflag as output. For now only indicates if data aren't equidistant within 1%
;   5aug2020   added init_minflux kw. The default for the initial min-time estimate is now the middle of input lc (better for noisy data)
;  17nov2020  replaced int() function calls with fix() 
;  18nov2020  uses now linfit() for linear transform between foldidf and time (more robust against scatter in time vs foldidf)
;  30nov2020  replaced isa() with keyword_set() statements. Added some debug/diagnostics stuff
;   4feb2021  added plotting of id numbers in lightcurve plot, and the noplnum kw
;  15nov2023  minor updates in program description in header, code not modified
;  21nov2023  some more printing if debug=2 or 3 ; updates in description
;  23nov2023  fixed the citation of Deeg 2021 paper


;Citing this code: Deeg, H.J. 2021, "A Modified Kwee-Van Woerden Method for Eclipse Minimum 
;                     Timing with Reliable Error Estimates", Galaxies, vol. 9, issue 1, p. 1 
  
;COPYRIGHT (C) 2020 Hans J. Deeg;   
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    any later version.
;
;    As a special exception, this program may be compiled with the
;    Interactive Data Language (IDL) and be linked to libaries
;    pertaining to IDL.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program in file 'licence.txt'; if not, write to
;    the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
;    Boston, MA 02110-1301 USA
;----------------------------------

  errflag=0
  if ~keyword_set(debug) then debug =0 ;debugflag
  if ~keyword_set(nfold) then nfold = 5 ; default number of foldings
  npts=n_elements(time)
  minid=npts/2                  ;ID of approximate center (and minimum) of lc

  if ~keyword_set(notimeoff) then begin ;offset times by value of central pt (lowered to preceeding increment of 0.1 )
     time0=floor(time[minid]*10D,/L64)/10D ;offset in time to first decimal          
     time-=time0
     tmpstr=string(time0,format='(F12.1)') 
     xtitlestr='time + '+tmpstr.trim()
  endif else begin              ;use original time values
     time0=0D
     xtitlestr='time'
  endelse

  if debug ge 3 then begin
     for i=0,npts-1 do begin
        print,i, flux[i]
     endfor
  endif
  if  ~keyword_set(noplot) then begin
     p1=plot(time,flux,'o:',xtitle=xtitlestr,ytitle='flux',title='Input lightcurve (Pts with crosses are used in KvW)')
     if ~keyword_set(noplnum) then begin ;plot ID-numbers of data points 
        for i=0,npts-1 do begin
           t1=text(time[i],flux[i],' '+strcompress(string(i),/rem),/data) 
        endfor
     endif
  endif
;test if data are equidistant within 1%
  difftime=time[1:*]-time[0:-2]
  mediandiff=median(difftime)
  if ( max(difftime) - mediandiff ) / mediandiff ge 0.01 or ( mediandiff-min(difftime)) / mediandiff ge 0.01 then begin
                                ;if ~keyword_set(noplot) then p=plot(time,flux,'+')
     print,'KVW: WARN: Time-points not equidistant within 1%. Raised errorflag '
     errflag +=1
  endif

  if keyword_set(init_minflux) then minflux=min(flux, minid) ;use as preliminary minid the flux minimum
  if debug ge 1 then print,'minid: ',minid



  noffr=(nfold-1)/4.            ;max offset in folding indices around minid 
  noff=fix(noffr)    


  minfoldid=minid-noff          ;min and max of IDs on which we fold (integer numbers)
  maxfoldid=minid+noff  
  nleft=minfoldid               ;need to test this with assymnetric lcs in either side
  nright=npts-maxfoldid-1  
  Z=nleft<nright                ;number of flux-pairs that are formed
  minpid=minfoldid-Z
  maxpid=maxfoldid+Z


  if debug ge 2 then print,'points considered in pairings:',minpid, ' to ',maxpid
  if  ~keyword_set(noplot) then p1=plot(time[minpid:maxpid],flux[minpid:maxpid],'+',/over)
  if debug ge 1 then  print,'Z: ' ,Z
  if Z lt 3 then message,'Error: Less than 3 points in in/egress can be paired. Decrease nfold parameter or provide more datapoints for in- or egress'


  S=fltarr(nfold)               ;sum of squares
  segid=0
  foldidf=fltarr(nfold)
  for segid=0,nfold-1 do begin               ;main loop performing foldings
     foldidf[segid]= minid-noffr+segid/2.    ;fractional foldid:  mindid-noff, mindid-noff+0.5,..minid,..mindid+noff-0.5, mindid+noff
     if debug ge 3 then print,'folding on ',foldidf[segid]
     foldidi = fix(foldidf[segid]+0.0001) ; integer value of foldid
     S[segid]=0.
     if abs (foldidf[segid] - foldidi )  le 0.01 then begin ; foldidf  is integer (folding on a point)         
        for i=1,Z do begin
           idlo=foldidi-i
           idhi=foldidi+i
           if debug ge 2 then print,'even',idlo,idhi
           S[segid] +=(flux[idlo]-flux[idhi])^2
        endfor
     endif else begin           ;foldidf is int.+half (folding between points)
        for i=1,Z do begin
           idlo=foldidi-i+1
           idhi=foldidi+i
           if debug ge 2 then print,'odd',idlo,idhi
           S[segid] +=(flux[idlo]-flux[idhi])^2
        endfor
     endelse
  endfor
  
;define a linear transformation between the values of foldidf and time: time = cl[0] + cl[1] * foldidf
  minidf=fix(foldidf[0])        ;get next integer below and above the range of foldidf (spans at least 3 integers)
  maxidf=fix(foldidf[-1]+0.50001)
  cl=linfit([minidf:maxidf],time[minidf:maxidf],/double )
  if debug ge 2 then print,'cl: ',cl
  timef=dblarr(nfold)           ;create time points at id-values of foldings
  for i=0,nfold-1 do timef[i]=cl[0]+cl[1]*(foldidf[i])
  if debug ge 2 then print,'timef (orig): ',timef

;test if timef is OK
;for i=0, nfold-1 do print,foldidf[i],timef[i]

  if  ~keyword_set(noplot) then p3=plot(timef,S,'o:',ytitle='S',xtitle=xtitlestr,title='KvW S-values of foldings and fit')       
  if  ~keyword_set(noplot) then  p3tax=axis('X',COORD_TRANSFORM=[-cl[0]/cl[1],1/cl[1]],tickinterval=0.5,minor=0,title='fold-ID') ;reverse transform to foldidf

;check if S has an imbalance relative to min(S) larger than one pt to either side, and cut it down to this
  if nfold ge 5 then begin
     minS=min(S,minSid)
     nSleft=minSid
     nSright=nfold-minSid-1
     if nSleft - nSright ge 2 then begin      ; >=2 points more on left than on right; cut left pts
        timef=timef[nSleft - nSright-1:* ]    ;timef: are time-values that are used in fit of S
        S=S[nSleft - nSright-1:* ]
     endif
     if nSright - nSleft ge 2 then begin ; >=2 points more on right than on left; cut right pts
        timef=timef[0:nSleft-nSright ]
        S=S[0:nSleft-nSright]
     endif
     if debug ge 2 then print,'nSleft: ',nSleft,' nSright: ',nSright

  endif                                                      ;nfold ge 5
  if  ~keyword_set(noplot) then  p3=plot(timef,S,'+:',/over)    ;plot pts used in polyfit with cross

    if debug ge 2 then print,'foldidf', foldidf
    if debug ge 2 then print,'timef', timef
    if debug ge 2 then print,'      S',S

  
  cp=poly_fit(timef,S,2,/double,status=status) ;fit 2nd-order polynome
  if debug ge 1 then print,'cp (kvw) = ',cp
;overplot the fitted parabole
  if  ~keyword_set(noplot) then begin
     nfitp=50                   ;number of pts to plot
     Sfit=fltarr(nfitp)
     idfit=fltarr(nfitp)
     for j=0,nfitp-1 do begin
        idfit[j]= timef[0] + (timef[-1]-timef[0])*j/(nfitp-1.) 
        Sfit[j]=cp[0]+cp[1]*idfit[j]+cp[2]*idfit[j]^2
     endfor
     p4=plot(idfit,Sfit,/over)
  endif                         ;if  ~keyword_set(noplot) 

  mintim=-cp[1]/(2.*cp[2])  +time0                             ;eq 3 from KvW56
  kvwminerr=sqrt((4*cp[2]*cp[0]-cp[1]^2)/(4*cp[2]^2*(Z-1)))    ;eq 4 from KvW56   ;original error calc
  if debug ge 3 then print,format='(a,f10.6)','cp[0] orig: ',cp[0]
  if debug ge 3 then print,format='(a,f15.7,a,f9.7)','mintim_orgKvW: ',mintim,'+-',kvwminerr
  if debug ge 3 then begin
     minS_org = cp[0]-cp[1]^2/(4*cp[2]) ;min value of fitted parabole; eq. 2 from KvW
     cp0org=cp[0]
     print,'minS_org:   ',minS_org
  endif


;determine flux-rms if not supplied by kw 
  if ~(keyword_set(rms)) then begin   
     rms=sqrt(min(S)/(2*(Z-1))) ;automatic rms determination from lowest S at any folding
     if ~keyword_set(noprint) then print,'KVW: flux-rms from best folding: ',rms
  endif

;revised calculation of timing error
  cp[0]=(Z-1)*2*rms^2 +cp[1]^2/(4*cp[2]) ;recalc cp[0] so that minS = (Z-1)*2*rms^2
  minerr=sqrt((4*cp[2]*cp[0]-cp[1]^2)/(4*cp[2]^2*(Z-1))) ;eq 4 from KvW56   ;using revised cp[0]
  if debug ge 1 then begin
     print,'cp (revised): ',cp
     print,'rms ',rms
     print,'minerr: ',minerr
  endif
  
  if ~keyword_set(noprint) then print,format='(a,f15.7,a,f9.7)','KVW: mintime: ',mintim,'+-',minerr

  if debug ge 3 then begin
     print,'minerr_alt: ',rms*sqrt(2/cp[2])                    ;alternative eq. for minerr (my eq 8)
     minS = cp[0]-cp[1]^2/(4*cp[2])                            ;min value of fitted parabole; eq. 2 from KvW
     print,'minS:   ',minS
     print,'cp[0] -cp[0]org or minS-minSorg',cp[0]-cp0org
     print,'min chi^2: ',minS / (2.*rms^2) ;min value in chi^2 
  endif
;stop
  return,mintim
end
