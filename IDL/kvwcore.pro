function kvwcore, time, flux, minerr=minerr, kvwminerr=kvwminerr,  init_minflux=init_minflux, rms=rms,  nfold=nfold, errflag=errflag
  
;kvwcor.pro is a plugin replacement for kvw.pro in which non-essential code (options for graphics, printing, debugging, time-offset)
;has been removed

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
                                ;errflag output flag with error-codes. 0: OK, 1: data not equidistant spaced
  
;HJD 19nov2019 initial version
;  24jan2020   introduced internal offset of input times by value of central time-point (floored to first digit) to avoid numerical
;               problems with polyfit on large time-values (even values close to 1 are problematic!). Also new kw notimeoff 
;  31jul2020   added errflag as output. For now only indicates if data aren't equidistant within 1%
;   5aug2020   added init_minflux kw. The default for the initial min-time estimate is now the middle of input lc (better for noisy data)
;  17nov2020  replaced int() function calls with fix() 
;  18nov2020  uses now linfit() for linear transform between foldidf and time (more robust against scatter in time vs foldidf)
;  18nov2020  kvwcore.pro generated from kvw.pro
;  30nov2020  replaced isa() with keyword_set() statements, similar to today's version of kvw.pro
;  21nov2023   updates in program description in header, code not modified except for order of kw's
;  23nov2023  fixed the citation of Deeg 2021 paper

;Citing this code: Deeg, H.J. 2021, "A Modified Kwee-Van Woerden Method for Eclipse Minimum 
;                     Timing with Reliable Error Estimates"Galaxies, vol. 9, issue 1, p. 1 
  
;COPYRIGHT (C) 2020, 2023 Hans J. Deeg;   
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
  if ~keyword_set(nfold) then nfold = 5 ; default number of foldings

  npts=n_elements(time)
  minid=npts/2                               ;ID of approximate center (and minimum) of lc
  time0=floor(time[minid]*10D,/L64)/10D      ;offset in time to first decimal          
  time-=time0

  
;test if data are equidistant within 1%
  difftime=time[1:*]-time[0:-2]
  mediandiff=median(difftime)
  if ( max(difftime) - mediandiff ) / mediandiff ge 0.01 or ( mediandiff-min(difftime)) / mediandiff ge 0.01 then begin
     print,'KVW: WARN: Time-points not equidistant within 1%. Raised errorflag '
     errflag +=1
  endif

  if keyword_set(init_minflux) then minflux=min(flux, minid) ;use as preliminary minid the flux minimum

  noffr=(nfold-1)/4.            ;max offset in folding indices around minid 
  noff=fix(noffr)    

  minfoldid=minid-noff          ;min and max of IDs on which we fold (integer numbers)
  maxfoldid=minid+noff  
  nleft=minfoldid               ;need to test this with assymnetric lcs in either side
  nright=npts-maxfoldid-1  
  Z=nleft<nright                ;number of flux-pairs that are formed
  minpid=minfoldid-Z
  maxpid=maxfoldid+Z


  if Z lt 3 then message,'Error: Less than 3 points in in/egress can be paired. Decrease nfold parameter or provide more datapoints for in- or egress'


  S=fltarr(nfold)               ;sum of squares
  segid=0
  foldidf=fltarr(nfold)
  for segid=0,nfold-1 do begin               ;main loop performing foldings
     foldidf[segid]= minid-noffr+segid/2.    ;fractional foldid:  mindid-noff, mindid-noff+0.5,..minid,..mindid+noff-0.5, mindid+noff
     foldidi = fix(foldidf[segid]+0.0001)    ; integer value of foldid
     S[segid]=0.
     if abs (foldidf[segid] - foldidi )  le 0.01 then begin ; foldidf  is integer (folding on a point)         
        for i=1,Z do begin
           idlo=foldidi-i
           idhi=foldidi+i
           S[segid] +=(flux[idlo]-flux[idhi])^2
        endfor
     endif else begin           ;foldidf is int.+half (folding between points)
        for i=1,Z do begin
           idlo=foldidi-i+1
           idhi=foldidi+i
           S[segid] +=(flux[idlo]-flux[idhi])^2
        endfor
     endelse
  endfor
  
;define a linear transformation between the values of foldidf and time: time = cl[0] + cl[1] * foldidf
  minidf=fix(foldidf[0])        ;get next integer below and above the range of foldidf (spans at least 3 integers)
  maxidf=fix(foldidf[-1]+0.50001)
  cl=linfit([minidf:maxidf],time[minidf:maxidf],/double )

  timef=dblarr(nfold)           ;create time points at id-values of foldings
  for i=0,nfold-1 do timef[i]=cl[0]+cl[1]*(foldidf[i])


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
  endif                         ;nfold ge 5
  
  cp=poly_fit(timef,S,2,/double,status=status) ;fit 2nd-order polynome

  mintim=-cp[1]/(2.*cp[2])  +time0                             ;eq 3 from KvW56
  kvwminerr=sqrt((4*cp[2]*cp[0]-cp[1]^2)/(4*cp[2]^2*(Z-1)))    ;eq 4 from KvW56   ;original error calc

;determine flux-rms if not supplied by kw 
  if ~(keyword_set(rms)) then begin   
     rms=sqrt(min(S)/(2*(Z-1))) ;automatic rms determination from lowest S at any folding
  endif

;revised calculation of timing error
  cp[0]=(Z-1)*2*rms^2 +cp[1]^2/(4*cp[2])                    ;recalc cp[0] so that minS = (Z-1)*2*rms 
  minerr=sqrt((4*cp[2]*cp[0]-cp[1]^2)/(4*cp[2]^2*(Z-1)))    ;eq 4 from KvW56
  
  return,mintim
end
