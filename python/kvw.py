import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def kvw(time, flux, init_minflux=False, rms=None, nfold=5, notimeoff=False,  noplot=False, noprint=False, debug=0):
    """
;Returns eclipse mid (minimum)-time using the Kwee Van-Woerden (1956) method with revised timing error, following Deeg 2020 (Galaxies, vol. 9, issue 1, p. 1)
;Data-points that are equidistant in time are required and the lightcurve should only contain the eclipse (no off-eclipse data).
;For the initial guess of the minimum time, two options are given, controlled by the keyword init_minflux : By default (init_minflux=0), the middle of 
;the lightcurve is assumed. Else (init_minflux=1), the time of the flux-minimum is used. This option is fine or preferential for low-noise lightcurves showing 
;a clear flux-mimimum not extending over more than 2-3 points. 
;If the rms (point-to-point noise) of the off-eclipse flux is known, it should be supplied by keyword rms. Else, rms is estimated from S of the best 
; flux-pairing, which is assumed to be dominated by flux measurement errors. For details on calculation of that rms, see Deeg 2020. 
;By default, the time of minimum is derived from 5 pairings (nfold paramter) that fold at i-1,i-0.5,i,i+0.5,i+1, where i is the index of the point of minimum 
; flux. The original KvW algorithm uses nfold=3.

;input 
    ;time   vector with time values in ascending order
    ;flux   vector with corresponding flux values
;output 
    ;tuple with 4 values: 
    ;time of mimimum, error of minimum (method by Deeg 2020), error of minimum (KvW's orignal method), and error-code.
    ;error-codes: 0: OK, data are equidistant within 1% against medium spacing.  1: data not equidistantly spaced

;keywords
    ;init_minflux input:  By default, the middle of the lightcurve is used as intial guess of the minimum time. If 1, uses as initial min-time 
    ;                     estimate the value of lowest flux.
    ;rms     input: average measurement error of individual flux values (in units of the flux values). rms should usually be supplied.     
    ;nfold   input: number of foldings on which to perform pairings of flux values, default=5
    ;notimeoff input: if set to 1, disables internal offsetting of time values. This is save to use if eclipse min-time is within +-0.1 of time of the intial min-time estimate)
    ;noplot  input: if set to 1, supresses plots
    ;noprint input: if set to 1, supresses all text-output
    ;debug   input: debugging flag. If 1, prints some, if 2, prints lots, if 3 even more intermediate values
 
  HJD 16nov2023:  First version of kvw.py, translated from kvw.pro version 15nov2023 using Chatgpt4. 
      21nov2023: Revised code that delivers identical numerical results as IDL code kvw.pro (ecxept for graphics, which are simpler)

;
;CITING this code: Deeg, H.J. 2020, "A Modified Kwee-Van Woerden Method for Eclipse Minimum 
;                     Timing with Reliable Error Estimates"Galaxies, vol. 9, issue 1, p. 1     
      
;COPYRIGHT (C) 2023 Hans J. Deeg;   
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
    """
    np.seterr(invalid='ignore') #supress warnings on negative sqrt
    errflag = 0
    npts = len(time)
    minid = npts // 2

    if not notimeoff:
        time0 = np.floor(time[minid] * 10) / 10
        time -= time0
        xtitlestr = f'time + {time0:.1f}'
    else:
        time0 = 0
        xtitlestr = 'time'

    if debug >= 3:
        for i in range(npts):
            print(i, flux[i])

    if not noplot:
        plt.figure()
        plt.plot(time, flux, 'o:',fillstyle='none')
        plt.xlabel(xtitlestr)
        plt.ylabel('flux')
        plt.title('Input lightcurve (Pts with crosses are used in KvW)')
        plt.show()

    # Check if data are equidistant within 1%
    difftime = time[1:] - time[:-1]
    mediandiff = np.median(difftime)
    if np.max(difftime) / mediandiff >= 1.01 or mediandiff / np.min(difftime) >= 1.01:
        print('KVW: WARN: Time-points not equidistant within 1%. Raised errorflag')
        errflag += 1

    if init_minflux:
        minflux = np.min(flux)
        minid = np.argmin(flux)
        if debug >=2:
            print('minid from init_minflux',minid)


    if debug >= 1:
        print('minid ',minid)
  
    noffr = (nfold - 1) / 4.0
    noff = int(noffr)
    minfoldid = minid - noff
    maxfoldid = minid + noff
    nleft = minfoldid
    nright = npts - maxfoldid - 1
    z = min(nleft, nright)
    minpid = minfoldid - z
    maxpid = maxfoldid + z

    if debug >= 2:
        print('points considered in pairings:',minpid, ' to ',maxpid)
    
    if not noplot:
        plt.plot(time[minpid:maxpid+1], flux[minpid:maxpid+1], 'x',fillstyle='none')
    if debug >= 1:
        print('Z: ',z)

    if z < 3:
        raise ValueError('Error: Less than 3 points in in/egress can be paired. Decrease nfold parameter or provide more datapoints for in- or egress')

    s = np.zeros(nfold)
    foldidf = np.zeros(nfold)
    for segid in range(nfold):
        foldidf[segid] = minid - noffr + segid / 2.0
        foldidi = int(foldidf[segid] + 0.0001)
        if abs(foldidf[segid] - foldidi) <= 0.01:
            for i in range(1, z + 1):
                idlo = foldidi - i
                idhi = foldidi + i
                s[segid] += (flux[idlo] - flux[idhi]) ** 2
        else:
            for i in range(1, z + 1):
                idlo = foldidi - i + 1
                idhi = foldidi + i
                s[segid] += (flux[idlo] - flux[idhi]) ** 2

    # Linear transformation between foldidf and time
    minidf = int(np.floor(foldidf[0]))
    maxidf = int(np.floor(foldidf[-1] + 0.50001))
    popt, _ = curve_fit(lambda x, a, b: a + b * x, np.arange(minidf, maxidf + 1), time[minidf:maxidf + 1])
    if debug >= 2:
        print('popt: ',popt)
    timef = popt[0] + popt[1] * foldidf
    if debug >=2:
        print('timef (orig): ',timef)

    if not noplot:
        plt.figure()
        plt.plot(timef, s, 'o:', fillstyle='none')
        plt.xlabel(xtitlestr)
        plt.ylabel('S')
        plt.title('KvW S-values of foldings and fit')
        plt.show()


# check if S has an imbalance relative to min(S) larger than one pt to either side, and cut it down to this
    if nfold >= 5:
        minS = np.min(s)
        minSid = np.argmin(s)
        nSleft = minSid
        nSright = nfold - minSid - 1

        if nSleft - nSright >= 2:  # More than 2 more points on left than on right, cut left points
            timef = timef[nSleft - nSright - 1:]
            s = s[nSleft - nSright - 1:]

        if nSright - nSleft >= 2:  # More than 2 more points on right than on left, cut right points
            timef = timef[:nSleft - nSright+1]
            s = s[:nSleft - nSright+1]
        if debug >= 2:
            print('nSleft: ',nSleft,' nSright: ',nSright)
    
    if not noplot:
        plt.plot(timef, s, 'x',fillstyle='none')

    if debug >= 2:
        print('foldidf ', foldidf)
        print('timef ',timef)
        print('     S ',s)

    # Fit a second-order polynomial
#    cp, _ = np.polyfit(timef, s, 2, full=False)
    cp = np.polyfit(timef, s, 2, full=False)
    if debug >= 1:
        print('cp (kvw): ',cp)

    if not noplot:  #plot fitted polynomial
        nfitp=50    #number of pts to plot
        Sfit= np.zeros(50, dtype=float)
        idfit= np.zeros(50, dtype=float)
        for j in range (0,nfitp):
            idfit[j]= timef[0] + (timef[-1] - timef[0]) * j / (nfitp - 1.) 
            Sfit[j]=cp[2] + cp[1] * idfit[j] + cp[0] * idfit[j] ** 2
        #print('idfit ',idfit)
        #print('Sfit',Sfit)  
        plt.plot(idfit,Sfit)          

    mintim = -cp[1] / (2.0 * cp[0]) + time0
    kvwminerr = np.sqrt((4 * cp[0] * cp[2] - cp[1] ** 2) / (4 * cp[0] ** 2 * (z - 1)))

    # Determine flux-rms if not supplied
    if rms is None:
        rms = np.sqrt(np.min(s) / (2 * (z - 1)))

    # Revised calculation of timing error
#    cp[0] = (z - 1) * 2 * rms ** 2 + cp[1] ** 2 / (4 * cp[2])
#    minerr = np.sqrt((4 * cp[2] * cp[0] - cp[1] ** 2) / (4 * cp[2] ** 2 * (z - 1)))
    cp[2] = (z - 1) * 2 * rms ** 2 + cp[1] ** 2 / (4 * cp[0])
    #cp[2]= 38.583897
    minerr = np.sqrt((4 * cp[0] * cp[2] - cp[1] ** 2) / (4 * cp[0] ** 2 * (z - 1)))
    if debug >= 1:
        print('cp (revised): ',cp)
        print('type of cp:', type(cp))
        print('rms:',rms)
        print('minerr: ', minerr)
    
    if not noprint:
        print(f"mintime: {mintim:15.7f} +- {minerr:9.7f}")
    return mintim, minerr, kvwminerr, errflag

# Example usage of the function
# time = np.array([...])  # Your time data
# flux = np.array([...])  # Your flux data
# result = kvw(time, flux)




# Function to read the data from the file
def read_lightcurve(filename):
    time, flux = [], []
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                parts = line.split()
                time.append(float(parts[0]))
                flux.append(float(parts[1]))
    return np.array(time), np.array(flux)


#
"""
; Demonstrator of KvW.pro to determine min-times of two eclipses of CM Dra
; in TESS light curves, with revised KvW method by Deeg 2020.
; The two light curves are:
; - Primary eclipse at epoch 7024, generating also Fig 1 and first entry
;    in Table 1  
; - Incomplete primary at epoch 7023, generating also Figs.3 and 4
; Printed output gives also the timing error from the orignal equation by 
; Kvw 1956  (with a NaN at epoch 7024).
;
;both curves were extracted from PDCSAP_FLUX of MAST file
;tess2019253231442-s0016-0000000199574208-0152-s_lc.fits
;and processed as described in the paper (Deeg 2020)

;Executing this code,  (in ipython: 'run kvw') 
;the text-output should be:
;infile= CMDra7024.lc
;mintime:   58739.9291169+-0.0000125 orig. KvW error:      -NaN
;----------------------------------
;infile= CMDra7023.lc
;mintime:   58738.6607358+-0.0000191 orig. KvW error: 0.0000662
;----------------------------------
"""
# Path to the file (you should adjust this to your files' location)
inpath='../example_data/'
infile = [inpath+'CMDra7024.lc',inpath+'CMDra7023.lc']
#infile = [inpath+'testleft.lc',inpath+'testirreg.lc']   #test-curves derived from CMDra7024.lc: For left-sided eclipse coverage and for irregular spacing

for lc in infile:

    # Read the data
    time, flux = read_lightcurve(lc)

    #rms of off-eclipse data, determined from revision of all Sector 16  CM Dra eclipses, see paper Deeg 2020
    rms=0.00138   
    
    # Apply the kvw function
    mintim, minerr, kvwminerr, errflag = kvw(time, flux, nfold=5, init_minflux=1, rms=rms, noplot=0, noprint=1,debug =0)
    
    # Output the results
    print('infile= ',lc )
    print(f"mintime: {mintim:15.7f} +- {minerr:9.7f} orig. KvW error: {kvwminerr:9.7f}")
    #print('errflag: ',errflag)
    print('----------------------------------')