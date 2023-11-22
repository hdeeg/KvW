pro rdtab, infilename,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,vtype=vtype,comc=com,headc=head,headarr=headstr,firstrow=firstrow,delstr=delstr,minlen=minlen,verbose=verbose,ndat=ndat,debugcol=debugcol,maxdatrow=maxdatrow, maxheadrow=maxheadrow, datarr=datarr, maxcol=maxcol, nowarn=nowarn
;flexible reader for alphanumeric tables extracting up to 30 columns
;reads as many columns as are passed by parameters c1 .. c30
;optionally reads in header-lines or ignores comment lines
;by default, space separated values are returned as float-arrays

;keyword vtype is array with numbers indicating variable type. Common
;            ones:  0 ignore column (no coresponding variable c1..c20
;            returned), 1 byte, 2 int, 3 long, 4 float, 5 double, 
;            7 string. For other ones, see IDL-help on  SIZE
;            if vtype is not given, double is used as default 
;keyword comc is a character indicating start of comment lines. Comment
;        lines are ignored and not returnd. If this
;        keyword is not set, all lines are considered data-lines
;keyword headc is character indicating start of a header lines. Header
;        lines are returned in array 'headarr'
;keyword headarr indicates a string array in which header-lines are returned
;keyword firstrow is number of first row to read, starting at 1. Default is 1. Lines previous to this are ignored completely, be they header or not.
;keyword maxdatrow: if not set, rdtab will read the entire table up to a max of 10^6 data lines. If set, rdtab will only read maxdatrow data lines (or to the end of the table, whatever comes first).
;keyword maxheadrow: if not set, rdtab will read the entire header up to max of 2000 lines. If set, rdtab will read up to maxheadrow data lines.
  
;keyword delstr is string that deliminites columns, default is  " ". See 'pattern' in help of strsplit
;keyword minlen gives min number of characters in valid data lines, default=1

;keyword verbose gives some feedback on number of lines red
;keyword ndat returns number of data-rows as output to a defined variable
;keyword datarr if set, returns data in ndat*maxcol -sized string array. 
;        maxcol is number of columns in first data-line or set by maxcol keyw. 
;        datarr must be a defined input variable of any type or value.
;        no values are assigned to variables c0...c30..
;keyword maxcol   only used together with datarr keyw, gives number of columns to be read into datarr. .
;keyword nowarn   supresses warning messages

;example 1: rdtab,'infile.txt',a,y,vtype=[7,0,2],delstr='&',comc='c',headc='#',headarr=headarr
;  reads from infile.txt the first column as string into variable 'a', ignores
;  the second column, and reads the third column as integer into 'y', using as
;  delimiter the '&' chracters as in a tex table. It ignores lines starting 
;  with 'c' and returns lines starting with '#' in an
;  array of header-lines called headarr. 
;
;example 2: 
;  outarr=0  ;outarr needs to be a defined variable
;  rdtab,'infile.csv',delstr=';',datarr=outarr,firstrow=2,maxcol=4
;         reads a semi-colon delimited csv file into a string array 'outarr', with 4 columns, 
;         and ignoring the first row.


;11nov03 HJD based on readgentab.pro
;18oct05 HJD insertion of vtype, headc, headarr keywords
;12sep07 HJD inserted ndat output keyword; doesn't crash anymore on
;files without data-rows
;20may08 HJD extended to 13 data columns
;27aug08 HJD extended to 20 data columns
;10sep08 HJD fixes to work with long tables (up to 1e6 data lines)
;12nov08 HJD fixed  bug incomplete reads on long tables (set nl,ncom,nhead variables to long)
;21nov08 HJD added minlen keyword
;21nov08 HJD allowed ignoring of input table colums by setting zero in vtype array
;11may10 HJD extended to 25 data columns. keyword debugcol prints column which is helpful to find type-conversin errors.
;29may10 HJD extended to 30 data columns.
;17sep13 HJD inserted keyw maxdatrow and maxheadrow to override the defaults of 10^6 data-lines or  2000 header lines that can maximally be ingested
;21oct2015 HJD Changed default data-type to double
;4mar2016 HJD redefined maxdatrow and maxheadrow to also mean the max number of rows that are read.
;23oct2016 HJD added datarr and maxcol keywords, to return original array of data-tokens
;11feb2018 HJD rewording of error message to make it clearer
;21nov2019 HJD datarr can now be smaller or larger than number of tokens in line; added kw nowarn; modified calling example given above
;26dec1019 HJD added example 2 on csv file

;max number of data red are set in this line:
  if not keyword_set(maxdatrow) then maxdatrow=1000000 ;max number of data-lines
  datstr=strarr(maxdatrow)
  if not keyword_set(maxheadrow)  then maxheadrow=2000 ;max number of header-lines
  headstr=strarr(maxheadrow)
  if ~isa(datarr) then begin ;standard procedure, defines ncol and nvalcols  
     ncol=n_params()-1
     if ncol gt 30 then print,'rdtab.pro: too many columns specified!'

     if NOT(keyword_set(vtype)) then  vtype=intarr(ncol)+5 ;set to double-precision for default
     dcol=where(vtype ne 0)                                ;array with data row numbers of valid data
     nvalcols=n_elements(dcol)                             ;number of valid columns specified in vtype
     if ncol gt nvalcols then begin                        ;this may happen only if vtype kw is used
        print,'Number of columns to read may not be larger than corresponding specifications in vtype-keyword!'
        stop
     endif 
  endif                        ; ~isa(datarr)

  if keyword_set(com) then cflg=1 else begin 
     cflg=0                     ;don't read comments
     com="dummy"
  endelse

  if keyword_set(head) then hflg=1 else begin 
     hflg=0                     ;don't read headers
     head="dummy"
  endelse
  if not(keyword_set(firstrow)) then firstrow=0 
  if ~isa(delstr) then delstr=" "
  if keyword_set(verbose) then verb=1 else verb=0
  if not(keyword_set(minlen)) then minlen=1 

  ndat=0L & ncom=0L &nhead=0L & nl=1L & tstr= ' '
  openr,UN1,/get_lun,infilename  
  while(not eof(UN1)) do begin
     readf,UN1,tstr
     if nl ge firstrow then begin
        case 1 of
           (cflg and (strmid(tstr,0,1) eq com)) or strlen(tstr) lt minlen: ncom=ncom+1 ;ignore comment lines or empty/short lines
           hflg and (strmid(tstr,0,1) eq head): begin
              if ndat lt maxheadrow then begin ;limit to maxheadrow header lines
                                ;if hstrip then tstr=strmid(tstr,1) ;strips the header-character
                 headstr(nhead)=tstr ;reads in a header-string
                 nhead=nhead+1
              endif
           end
           else: begin                        ;data are red
              if ndat lt maxdatrow then begin ;limit to maxdatrow data lines
                 datstr(ndat)=tstr 
                 ndat=ndat+1
              endif
           endelse
        endcase
     endif                      ;nl
     nl = nl+1
  endwhile

  free_lun,UN1
  if nhead ge 1 then headstr=headstr(0:nhead-1)
  if verb then print,ndat,' data lines red'
  if verb then print,ncom,' comment lines red'
  if verb then print,nhead,' header lines red'

  if delstr eq " " then presnull=0 else presnull=1 ; this sets an option for strsplit, 
                                ;so it defaults to splitting on spans of whitespace

  if ndat ge 1 then begin 
     datstr=datstr(0:ndat-1)
     if isa(datarr) then begin ;generate datarr output string array
        for k=0L,ndat-1 do begin
           tok=strsplit(datstr[k],delstr,/extract,preserve_null=presnull)
           ntok=n_elements(tok)
           if k eq 0 then begin 
              if keyword_set(maxcol) then noutcol=maxcol else noutcol =ntok ;number of columns in output arr to keyword value or to first data-line
              datarr=strarr(ndat,noutcol)
           endif    
           if ntok gt noutcol then if not(keyword_set(maxcol)) then $
              if ~keyword_set(nowarn) then tprint,'RDTAB: WARN: ',ntok,' tokens in data-line ',k,' versus ',noutcol,' tokens in data-array';more token in a data-line
           
           if ntok ge noutcol then datarr[k,*]=tok[0:noutcol-1] else begin 
               if ~keyword_set(nowarn) then print,'RDTAB: WARN:',ntok,' tokens in data-line ',k,' versus ',noutcol,' tokens in data-array' ;less token in a data-line
              datarr[k,0:ntok-1]=tok[0:ntok-1]
           endelse
 
        endfor                  ;k
     endif else begin           ;assign values to c1...c30 output arrays
;define data arrays
        c1=make_array(ndat,type=vtype[dcol[0]]) & if ncol le 1 then goto,enddef
        c2=make_array(ndat,type=vtype[dcol[1]])& if ncol le 2 then goto,enddef
        c3=make_array(ndat,type=vtype[dcol[2]])& if ncol le 3 then goto,enddef
        c4=make_array(ndat,type=vtype[dcol[3]])& if ncol le 4 then goto,enddef
        c5=make_array(ndat,type=vtype[dcol[4]])& if ncol le 5 then goto,enddef
        c6=make_array(ndat,type=vtype[dcol[5]])& if ncol le 6 then goto,enddef
        c7=make_array(ndat,type=vtype[dcol[6]])& if ncol le 7 then goto,enddef
        c8=make_array(ndat,type=vtype[dcol[7]])& if ncol le 8 then goto,enddef
        c9=make_array(ndat,type=vtype[dcol[8]])& if ncol le 9 then goto,enddef
        c10=make_array(ndat,type=vtype[dcol[9]])& if ncol le 10 then goto,enddef
        c11=make_array(ndat,type=vtype[dcol[10]])& if ncol le 11 then goto,enddef
        c12=make_array(ndat,type=vtype[dcol[11]])& if ncol le 12 then goto,enddef
        c13=make_array(ndat,type=vtype[dcol[12]])& if ncol le 13 then goto,enddef
        c14=make_array(ndat,type=vtype[dcol[13]])& if ncol le 14 then goto,enddef
        c15=make_array(ndat,type=vtype[dcol[14]])& if ncol le 15 then goto,enddef
        c16=make_array(ndat,type=vtype[dcol[15]])& if ncol le 16 then goto,enddef
        c17=make_array(ndat,type=vtype[dcol[16]])& if ncol le 17 then goto,enddef
        c18=make_array(ndat,type=vtype[dcol[17]])& if ncol le 18 then goto,enddef
        c19=make_array(ndat,type=vtype[dcol[18]])& if ncol le 19 then goto,enddef
        c20=make_array(ndat,type=vtype[dcol[19]])& if ncol le 20 then goto,enddef
        c21=make_array(ndat,type=vtype[dcol[20]])& if ncol le 21 then goto,enddef
        c22=make_array(ndat,type=vtype[dcol[21]])& if ncol le 22 then goto,enddef
        c23=make_array(ndat,type=vtype[dcol[22]])& if ncol le 23 then goto,enddef
        c24=make_array(ndat,type=vtype[dcol[23]])& if ncol le 24 then goto,enddef
        c25=make_array(ndat,type=vtype[dcol[24]])& if ncol le 25 then goto,enddef
        c26=make_array(ndat,type=vtype[dcol[25]])& if ncol le 26 then goto,enddef
        c27=make_array(ndat,type=vtype[dcol[26]])& if ncol le 27 then goto,enddef
        c28=make_array(ndat,type=vtype[dcol[27]])& if ncol le 28 then goto,enddef
        c29=make_array(ndat,type=vtype[dcol[28]])& if ncol le 29 then goto,enddef
        c30=make_array(ndat,type=vtype[dcol[29]])
        enddef: tmp=99  
;read in the main data block
        for k=0L,ndat-1 do begin
           tok=str_sep(strtrim(strcompress(datstr[k]),2),delstr)
           ntok=n_elements(tok)<ncol ;read all, or maximum of ncol columns

;finding bad entries --
           if keyword_set(debugcol) then begin ; print column to find type conv errors
              badcol=debugcol 
              print,k," ",tok[dcol[badcol]] 
           endif

           c1[k]=tok[dcol[0]]&if ntok le 1 then goto,endasig
           c2[k]=tok[dcol[1]]&if ntok le 2 then goto,endasig
           c3[k]=tok[dcol[2]]&if ntok le 3 then goto,endasig
           c4[k]=tok[dcol[3]]&if ntok le 4 then goto,endasig
           c5[k]=tok[dcol[4]]&if ntok le 5 then goto,endasig
           c6[k]=tok[dcol[5]]&if ntok le 6 then goto,endasig
           c7[k]=tok[dcol[6]]&if ntok le 7 then goto,endasig
           c8[k]=tok[dcol[7]]&if ntok le 8 then goto,endasig
           c9[k]=tok[dcol[8]]&if ntok le 9 then goto,endasig
           c10[k]=tok[dcol[9]]&if ntok le 10 then goto,endasig
           c11[k]=tok[dcol[10]]&if ntok le 11 then goto,endasig
           c12[k]=tok[dcol[11]]&if ntok le 12 then goto,endasig
           c13[k]=tok[dcol[12]]&if ntok le 13 then goto,endasig
           c14[k]=tok[dcol[13]]&if ntok le 14 then goto,endasig
           c15[k]=tok[dcol[14]]&if ntok le 15 then goto,endasig
           c16[k]=tok[dcol[15]]&if ntok le 16 then goto,endasig
           c17[k]=tok[dcol[16]]&if ntok le 17 then goto,endasig
           c18[k]=tok[dcol[17]]&if ntok le 18 then goto,endasig
           c19[k]=tok[dcol[18]]&if ntok le 19 then goto,endasig
           c20[k]=tok[dcol[19]]&if ntok le 20 then goto,endasig
           c21[k]=tok[dcol[20]]&if ntok le 21 then goto,endasig
           c22[k]=tok[dcol[21]]&if ntok le 22 then goto,endasig
           c23[k]=tok[dcol[22]]&if ntok le 23 then goto,endasig
           c24[k]=tok[dcol[23]]&if ntok le 24 then goto,endasig
           c25[k]=tok[dcol[24]]&if ntok le 25 then goto,endasig
           c26[k]=tok[dcol[25]]&if ntok le 26 then goto,endasig
           c27[k]=tok[dcol[26]]&if ntok le 27 then goto,endasig
           c28[k]=tok[dcol[27]]&if ntok le 28 then goto,endasig
           c29[k]=tok[dcol[28]]&if ntok le 29 then goto,endasig
           c30[k]=tok[dcol[29]]
           endasig: tmp=99 
        endfor                     
     endelse                    ; keyword_set(datarr)
  endif                         ;if ndat ge 1
                                ; end of reading stuff
;stop

end
