# KvW
Kwee Van Woerden (KvW) method for eclipse or transit minimum timing, with improved error calculation, as described in [Deeg (2021)](https://ui.adsabs.harvard.edu/abs/2020Galax...9....1D/abstract). 

The main programs are `kvw.py` (python) and `kvw.pro` (IDL).
The python and IDL versions deliver identical numerical results and very similar graphics output.

 
The code requires an input light curve of equidistant points that contains only the eclipse, without any off-eclipse points. A value for the rms (noise) of the input light curve is also requested (but not necessary). The eclipse minimum time is obtained using KvW's original method [(Kwee & Van Woerden 1956)](https://ui.adsabs.harvard.edu/abs/1956BAN....12..327K/abstract), but using more than the 3 reflections of KvW's original algorithm, with a default of 5 reflections (`nfold` parameter). The error of the minimum time is calculated following Deeg (2021); the error from KvW's original formula is also provided. Both IDL and phython codes are functions that return the eclipse minimum time with its error; they also provide optional text output, graphics, as well as several levels of debug information,.



## Execution of the demos:

python:
In directory with the codes, within IPython or similar interactive environment: `run kvw`

IDL:
From IDL cmd-line in directory with the codes:
`IDL> .rnew kvwdemo1`

The text-output should be in either language:
```
	infile= CMDra7024.lc
	mintime:   58739.9291169+-0.0000125 orig. KvW error:       NaN
	----------------------------------
	infile= CMDra7023.lc
	mintime:   58738.6607358+-0.0000191 orig. KvW error: 0.0000662
	----------------------------------
```
From `CMDra7024.lc,` the demos generate also Fig. 1 and the first entry in Table 1 of the paper. From `CMDra7023.lc`, the demos will generate Figs. 3 and 4. (Only the IDL version will generate exact reproductions)


## Package Content: 

`python` Directory with python code
- 'kvw.py'  code with the main kvw function. Executing `kvw.py` from command line will also provide a demo run, giving minimum times of the CM Dra lightcurves in example_data and generating some figures similar to the paper by Deeg 2021

`IDL`  Directory with IDL code 
- `kvw.pro` is the self-contained kvw code
- `kvwcore.pro` is a plugin replacement for `kvw.pro` in which non-essential code (options for graphics, printing, debugging, time-offset) has been removed. It is provided in order to facilitate translation into other languages.
- `kvwdemo1.pro` runs a demo of `kvw.pro`, providing minimum times of the lightcurves in example_data and generating some figures from the paper by Deeg 2021
- `rdtab.pro`  is a table-reader that used by kvwdemo1.pro


`example_data`:
- `CMDra7023.lc` and `CMDra7024.lc`  Lightcurves of CM Dra eclipses,
 	used by demos in kvw.py or kvwdemo1.pro.
	These curves are:
	- Primary eclipse at epoch 7024. 
	- Incomplete primary at epoch 7023;
	
	Both lightcurves were extracted from the PDCSAP_FLUX of the file `tess2019253231442-s0016-0000000199574208-0152-s_lc.fits`, available on NASA's MAST and processed as described in the paper (Deeg 2021).

- `*.lc`  Further lightcurves used in tests during development.
	
 
## Citing the KvW code
The preferred way is by citing [Deeg, H.J. 2021, "A Modified Kwee-Van Woerden Method for Eclipse Minimum Timing with Reliable Error Estimates", Galaxies, vol. 9, issue 1, p. 1](https://ui.adsabs.harvard.edu/abs/2020Galax...9....1D/abstract). 


