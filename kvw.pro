function kvw, time, flux, minerr=minerr, kvwminerr=kvwminerr, nfold=nfold, rms=rms, notimeoff=notimeoff,  init_minflux=init_minflux, noplot=noplot, noprint=noprint, errflag=errflag,debug=debug
  
;Returns eclipse mid (minimum)-time using the Kwee Van-Woerden (1956) method with revised timing error calc by Deeg+2020
;By default, the min-time is derived from 5 pairings that fold at i-1,i-0.5,i,i+0.5,i+1, where i is the index of the point of minimum 
; flux. This number can be modified with the nfold keyword (The origninal KvW algorithms uses nfold=3)
;If flux rms is known, it should be supplied by keyword rms. Else, rms is estimated from S of best flux-pairing, which is assumed to 
;be dominated by flux measurement errors. For details on error calculation, see Deeg+ 2020. 
;Data-points that are equidistant in time are required. 
  
;input 
      ;time   vector with time values in ascending order
      ;flux   vector with corresponding flux values
;output 
      ;  time of mimimum

;keywords
      ;minerr  output variable with error of the minimum time, method by Deeg+ 2020
      ;kvwminerr  output variable with error of the minimum time, original calculation by KvW 1956
      ;nfold   number of foldings on which to perform pairings of flux values, default=5
      ;rms     average measurement error of individual flux values (in units of the flux values)
      ;notimeoff  disables internal offsetting of time values. Is save to use if eclipse min-time is within +-0.1
      ;init_minflux  if set, uses as initial min-time estimate the value of lowest flux (instead of 
      ;the middle of the input lightcurve)
     ;noplot  supresses plots
      ;noprint supresses all printing
      ;errflag  flag with error-codes. 0: OK, bit 1: data not equidistant spaced
      ;debug   debugging flag that causes  printing of lots of stuff 
  
;HJD 19nov2019 initial version
;  24jan2020   introduced internal offset of input times by value of central time-point (floored to first digit) to avoid numerical
;               problems with polyfit on large time-values (even values close to 1 are problematic!). Also new kw notimeoff 
;  31jul2020   added errflag as output. For now only indicates if data aren't equidistant within 1%
;   5aug2020   added kw init_minflux. The default for the initial min-time estimate is now the middle of input lc (better for noisy data)



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
  if ~isa(nfold) then nfold = 5 ; number of foldings on which to perform pairings
  npts=n_elements(time)
  minid=npts/2   ;ID of approximate center (and minimum) of lc
  if ~keyword_set(notimeoff) then begin ;offset times by value of central pt (lowered to preceeding increment of 0.1 )
     time0=floor(time[minid]*10D,/L64)/10D   ;offset in time to first decimal          
     time-=time0
     tmpstr=string(time0,format='(F12.1)') 
     xtitlestr='time + '+tmpstr.trim()
  endif else begin              ;use original time values
     time0=0D
     xtitlestr='time'
  endelse

  if ~isa(debug) then debug =0  ;debugflag
  if debug then begin
     for i=0,npts-1 do begin
        print,i, flux[i]
     endfor
  endif
  if  ~keyword_set(noplot) then p1=plot(time,flux,'o:',xtitle=xtitlestr,ytitle='flux',title='KvW input lightcurve, crossed pts used')

;test if data are equispaced within 1%
  difftime=time[1:*]-time[0:-2]
  mediandiff=median(difftime)
  if ( max(difftime) - mediandiff ) / mediandiff ge 0.01 or ( mediandiff-min(difftime)) / mediandiff ge 0.01 then begin
     ;if ~keyword_set(noplot) then p=plot(time,flux,'+')
     print,'KVW: WARN: Time-points not equidistant within 1%. Raised errorflag '
     errflag +=1
endif

if keyword_set(init_minflux) then minflux=min(flux, minid)   ;use as preliminary minid the flux minimum



  noffr=(nfold-1)/4.            ;max offset in folding indices around minid 
  noff=int(noffr)    


  minfoldid=minid-noff          ;min and max of IDs on which we fold (integer numbers)
  maxfoldid=minid+noff  
  nleft=minfoldid               ;need to test this with assymnetric lcs in either side
  nright=npts-maxfoldid-1  
  Z=nleft<nright                ;number of flux-pairs that are formed
  minpid=minfoldid-Z
  maxpid=maxfoldid+Z


  if debug then print,'points considered in pairings:',minpid, ' to ',maxpid
  if  ~keyword_set(noplot) then p1=plot(time[minpid:maxpid],flux[minpid:maxpid],'+',/over)
  if debug then  print,'Z: ' ,Z
  if Z lt 3 then message,'Error: Less than 3 points in in/egress can be paired. Decrease nfold parameter or provide more datapoints for in- or egress'


  S=fltarr(nfold)               ;sum of squares
  segid=0
  foldidf=fltarr(nfold)
  for segid=0,nfold-1 do begin               ;main loop performing foldings
     foldidf[segid]= minid-noffr+segid/2.    ;fractional foldid:  mindid-noff, mindid-noff+0.5,..minid,..mindid+noff-0.5, mindid+noff
     if debug then print,'folding on ',foldidf[segid]
     foldidi = int(foldidf[segid]+0.0001) ; integer value of foldid
     S[segid]=0.
     if abs (foldidf[segid] - foldidi )  le 0.01 then begin ; foldidf  is integer (folding on a point)         
        for i=1,Z do begin
           idlo=foldidi-i
           idhi=foldidi+i
           if debug then print,'even',idlo,idhi
           S[segid] +=(flux[idlo]-flux[idhi])^2
        endfor
     endif else begin           ;foldidf is int.+half (folding between points)
        for i=1,Z do begin
           idlo=foldidi-i+1
           idhi=foldidi+i
           if debug then print,'odd',idlo,idhi
           S[segid] +=(flux[idlo]-flux[idhi])^2
        endfor
     endelse
  endfor
  if debug then print,'foldidf', foldidf
  if debug then print,'      s',S
  
;define a linear transformation between the values of foldidf and time: time = atr + btr * foldidf
;this could also be done between time and its own indices, but is prone to errors if there is some missing data-point 
  minidf=int(foldidf[0])        ;get next integer below and above the range of foldidf
  maxidf=int(foldidf[-1]+0.50001)
  ntimef=maxidf-minidf          ;number of pts of timef
  btr = (time[maxidf]-time[minidf])/ntimef
  atr = time[minidf]- btr*minidf

;test if transform is OK
;for i=0, ntimef-1 do print,i+minidf,time[i+minidf],atr+btr*(i+minidf)

  timef=dblarr(nfold)           ;create time points at id-values of foldings
  for i=0,nfold-1 do timef[i]=atr+btr*(foldidf[i])

;test if timef is OK
;for i=0, nfold-1 do print,foldidf[i],timef[i]

  if  ~keyword_set(noplot) then p3=plot(timef,S,'o:',ytitle='S',xtitle=xtitlestr,title='KvW S-values of foldings and fit')       
  if  ~keyword_set(noplot) then  p3tax=axis('X',COORD_TRANSFORM=[-atr/btr,1/btr],tickinterval=0.5,minor=0,title='fold-ID') ;reverse transform to foldidf

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
  endif                                                ;nfold ge 5
  if  ~keyword_set(noplot) then  p3=plot(timef,S,'+:',/over)    ;plot pts used in polyfit with cross
  
  cof=poly_fit(timef,S,2,/double,status=status) ;fit 2nd-order polynome

;overplot the fitted parabole
  if  ~keyword_set(noplot) then begin
     nfitp=50                   ;number of pts to plot
     Sfit=fltarr(nfitp)
     idfit=fltarr(nfitp)
     i=0
     for j=0,nfitp-1 do begin
        idfit[j]= timef[0] + (timef[-1]-timef[0])*j/(nfitp-1.) 
        Sfit[i]=cof[0]+cof[1]*idfit[j]+cof[2]*idfit[j]^2
        i+=1
     endfor
     p4=plot(idfit,Sfit,/over)
  endif                         ;if  ~keyword_set(noplot) 

  mintim=-cof[1]/(2.*cof[2])  +time0                               ;eq 3 from KvW56
  kvwminerr=sqrt((4*cof[2]*cof[0]-cof[1]^2)/(4*cof[2]^2*(Z-1)))    ;eq 4 from KvW56   ;original error calc
  if debug then print,format='(a,f14.6,a,f8.6)','mintim_org: ',mintim,'+-',kvwminerr

;determine flux-rms if not supplied by kw 
  if ~(keyword_set(rms)) then begin   
     rms=sqrt(min(S)/(2*(Z-1))) ;automatic rms determination from lowest S at any folding
     if ~isa(noprint) then print,'KVW: flux-rms from best folding: ',rms
  endif

;revised calculation of timing error
  cof[0]=(Z-1)*2*rms^2 +cof[1]^2/(4*cof[2])                     ;recalc cof[0] so that minS = (Z-1)*2*rms 
  minerr=sqrt((4*cof[2]*cof[0]-cof[1]^2)/(4*cof[2]^2*(Z-1)))    ;eq 4 from KvW56
  
  if ~isa(noprint) then print,format='(a,f14.6,a,f8.6)','KVW: mintime: ',mintim,'+-',minerr

  if debug then begin
     print,'minerr_alt: ',rms*sqrt(2/cof[2])                      ;alternative eq. for minerr (my eq 8)
     minS = cof[0]-cof[1]^2/(4*cof[2])                            ;min value of fitted parabole; eq. 2 from KvW
     print,'min S:   ',minS
     print,'min chi^2: ',minS / (2.*rms^2) ;min value in chi^2 
  endif
  return,mintim
end
