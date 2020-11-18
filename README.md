# KvW
Kwee Van Woerden (KvW) method for eclipse or transit minimum timing, with improved error calculation, as described in Deeg 2020 (submitted). The main program is kvw.pro, written in IDL. It requires an input light curve of equidistant points that contains only the eclipse, without any off-eclipse points. A value for the rms (noise) of the input light curve is also requested (but not necessary). The eclipse minimum time is obtained using KvW's original method, but the use of more than the 3 reflections of the original version is supported. The error of the minimum time is calculated following the formula in Deeg (2020), with an optional calculation using KvW's original formula. The principal output is the eclipse minimum time with its error, with different levels of auxiliary text and graphics output available.

Content: 
- kvw.pro is the self-contained kvw code
- kvwcore.pro is a plugin replacement for kvw.pro in which non-essential code (options for graphics, printing, debugging, time-offset) has been removed. It is provided in order to facilitate translation into other languages.
- 'Example' contains a demo wrapper to kvw.pro that generates some of the graphics and numerical output shown in the Deeg (2020) paper.  

