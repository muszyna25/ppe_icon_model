pro scores, plotdir, expnum, expref, inidate
;------------------------------------------------------------
; Plot scores from ICON experiments.
;
; run as: scores, '/fe1-daten/mkoehler/plots/icon', '82', '81', '20110101'
;
; Martin Koehler, Sep 2012
;------------------------------------------------------------

;print, 'Arguments: ', plotdir, ' ', expnum, ' ', expref, ' ', inidate
;print, ''

nlev=90
nlevpl=25
nlevzl=25
nword=18

close,/all
file1 = plotdir+"/nwp.exp"+expnum+"/scores_exp"+expnum+".txt"
file2 = plotdir+"/nwp.exp"+expref+"/scores_exp"+expref+".txt"
openr, 1, file1
openr, 2, file2
spawn,'wc -l '+file1, nlines1
spawn,'wc -l '+file2, nlines2
nlines1=strsplit(nlines1,/extract)
nlines2=strsplit(nlines2,/extract)
var1=make_array(nword,nlines1(0),/string)
var2=make_array(nword,nlines2(0),/string)
tmp=''
for nn=1,nlines1(0) do begin
  readf, 1, tmp
  var1(*,nn-1) = strsplit(tmp,/extract)
end
for nn=1,nlines2(0) do begin
  readf, 2, tmp
  var2(*,nn-1) = strsplit(tmp,/extract)
end

var2d  = ['TOT_PREC', 'TCC', 'TQV', 'TQC', 'TQI', 'ACCLHFL_S', 'ACCSHFL_S', $
         'ACCSOB_S', 'ACCTHB_S', 'ACCSOB_T', 'ACCTHB_T']
nvar2d = n_elements(var2d)

var3d  = ['T', 'U', 'V', 'P', 'QV', 'CC', 'CLWC', 'CIWC']   ;, 'CLWC', 'CRWC','CSWC'];'QTVAR', 'O3','P''QV', 'QC', 'QI'
nvar3d = n_elements(var3d)

var3pl = ['T', 'U', 'V', 'FI']   ;, 'Z'
nvar3pl= n_elements(var3pl)

var3zl = ['T', 'U', 'V', 'P', 'QC', 'QI']   ;, 'Z'
nvar3zl= n_elements(var3zl)

paperopenl
loadct,13

;--------------------------------------------------------------------------------------------

for nt=1,3 do begin

  CASE nt OF
    1: BEGIN
      step='24'
      nfor='10'
    END
    2: BEGIN
      step='120'
      nfor='6'
    END 
    3: BEGIN
      step='24'
      nfor='1'
    END 
  ENDCASE

;--------------------------------------------------------------------------------------------

  levtype='ml'

  xyouts, 0.15, 0.97, /normal, charsize=1.0, $
    'Comparison of exp'+expnum+' to exp'+expref+ $
    '   forecast hours: '+step+'   number of forecasts: '+nfor+ ' from '+inidate+'   '+levtype

  for nn=1,nvar3d do begin
   for nk=1,nlev do begin
   ;print,''
   ;print,'searching for variable ', var3d(nn-1)
    ind1 = where(var1(0,*) eq var3d(nn-1) and var1(2-1,*) eq levtype and var1(8-1,*) eq step and var1(10-1,*) eq nfor and var1(4-1,*) eq nk)
    ind2 = where(var2(0,*) eq var3d(nn-1) and var2(2-1,*) eq levtype and var2(8-1,*) eq step and var2(10-1,*) eq nfor and var2(4-1,*) eq nk)
   ;print, var3d(nn-1), ' L', nk, ' exp',expref, ' bias:', var2(16-1,ind2), 'rms:', var2(18-1,ind2)
    if ind1 eq -1 then print,'data1 doesnt exist: var ', var3d(nn-1), ' levtype ',levtype, ' step ', step, ' nfor ', nfor, ' nk ', nk
    if ind2 eq -1 then print,'data2 doesnt exist: var ', var3d(nn-1), ' levtype ',levtype, ' step ', step, ' nfor ', nfor, ' nk ', nk
    bias = abs(float(var1(16-1,ind1))) / ( abs(float(var1(16-1,ind1))) + abs(float(var2(16-1,ind2))) )
    rms  =     float(var1(18-1,ind1))  / (     float(var1(18-1,ind1))  +     float(var2(18-1,ind2))  )
    if ( abs(float(var1(16-1,ind1))) + abs(float(var2(16-1,ind2))) ) lt 1e-8 then begin
      bias = 0.5
    endif
    if (     float(var1(18-1,ind1))  +     float(var2(18-1,ind2))  ) lt 1e-8 then begin
      rms  = 0.5
    endif
  
    xxx =-0.11+nn*0.125
    yyy =0.93
    vert=0.01
  
    xyouts, xxx+0.03, yyy+0.015, var3d(nn-1), /normal, charsize=0.7
    header = 'level    rms      new'
    xyouts, xxx, yyy, header, /normal, charsize=0.7
    err1=var1(18-1,ind1)
    err2=var2(18-1,ind2)
    err =rms
    xoff=0.0
    amp =40.0
    scale = (err-0.5) * (err gt 0.5)
    xr= [0.0,0.01,0.01,0.0 ,0.0] * amp * scale(0) + xxx + xoff
    yr= [0.0,0.0 ,0.01,0.01,0.0] * 0.8            + yyy
    polyfill, xr   +0.04, yr   -nk*vert,              /normal, color=240
    scale=(0.5-err) * (err lt 0.5)
    xr=-[0.0,0.01,0.01,0.0 ,0.0] * amp * scale(0) + xxx + xoff
    polyfill, xr     +0.04, yr     -nk*vert,           /normal, color=160
    plots,    xr(3:4)+0.04, [yr(3)+0.002,yr(4)-0.002]-nk*vert,           /normal

    xyouts,   xr(0)     , yr(0)-nk*vert, strcompress(string(nk)), /normal, charsize=0.5
    xyouts,   xr(0)+0.07, yr(0)-nk*vert, err1,  /normal, charsize=0.4
   end
  end
  
  erase
  
;--------------------------------------------------------------------------------------------

  levtype='pl'

  xyouts, 0.15, 0.97, /normal, charsize=1.0, $
    'Comparison of exp'+expnum+' to exp'+expref+ $
    '   forecast hours: '+step+'   number of forecasts: '+nfor+ ' from '+inidate+'   '+levtype

  for nn=1,nvar3pl do begin
   for nk=1,nlevpl do begin
   ;print,''
   ;print,'searching for variable ', var3pl(nn-1)
    ind1 = where(var1(0,*) eq var3pl(nn-1) and var1(2-1,*) eq levtype and var1(8-1,*) eq step and var1(10-1,*) eq nfor and var1(4-1,*) eq nk)
    ind2 = where(var2(0,*) eq var3pl(nn-1) and var2(2-1,*) eq levtype and var2(8-1,*) eq step and var2(10-1,*) eq nfor and var2(4-1,*) eq nk)
   ;print, var3pl(nn-1), ' L', nk, ' exp',expref, ' bias:', var2(16-1,ind2), 'rms:', var2(18-1,ind2)
    bias = abs(float(var1(16-1,ind1))) / ( abs(float(var1(16-1,ind1))) + abs(float(var2(16-1,ind2))) )
    rms  =     float(var1(18-1,ind1))  / (     float(var1(18-1,ind1))  +     float(var2(18-1,ind2))  )
    if ( abs(float(var1(16-1,ind1))) + abs(float(var2(16-1,ind2))) ) lt 1e-8 then begin
      bias = 0.5
    endif
    if (     float(var1(18-1,ind1))  +     float(var2(18-1,ind2))  ) lt 1e-8 then begin
      rms  = 0.5
    endif
  
    xxx =-0.11+nn*0.14
    yyy =0.93
    vert=0.01
  
    xyouts, xxx+0.03, yyy+0.015, var3pl(nn-1), /normal, charsize=0.7
    header = 'level    rms      new'
    xyouts, xxx, yyy, header, /normal, charsize=0.7
    err1=var1(18-1,ind1)
    err2=var2(18-1,ind2)
    err =rms
    xoff=0.0
    amp =40.0
    scale = (err-0.5) * (err gt 0.5)
    xr= [0.0,0.01,0.01,0.0 ,0.0] * amp * scale(0) + xxx + xoff
    yr= [0.0,0.0 ,0.01,0.01,0.0] * 0.8            + yyy
    polyfill, xr   +0.04, yr   -nk*vert,              /normal, color=240
    scale=(0.5-err) * (err lt 0.5)
    xr=-[0.0,0.01,0.01,0.0 ,0.0] * amp * scale(0) + xxx + xoff
    polyfill, xr     +0.04, yr     -nk*vert,           /normal, color=160
    plots,    xr(3:4)+0.04, [yr(3)+0.002,yr(4)-0.002]-nk*vert,           /normal

    xyouts,   xr(0)     , yr(0)-nk*vert, strcompress(string(nk)), /normal, charsize=0.5
    xyouts,   xr(0)+0.07, yr(0)-nk*vert, err1,  /normal, charsize=0.4
   end
  end
  
  erase
  
;--------------------------------------------------------------------------------------------

  levtype='zl'

  xyouts, 0.15, 0.97, /normal, charsize=1.0, $
    'Comparison of exp'+expnum+' to exp'+expref+ $
    '   forecast hours: '+step+'   number of forecasts: '+nfor+ ' from '+inidate+'   '+levtype

  for nn=1,nvar3zl do begin
   for nk=1,nlevzl do begin
   ;print,''
   ;print,'searching for variable ', var3zl(nn-1)
    ind1 = where(var1(0,*) eq var3zl(nn-1) and var1(2-1,*) eq levtype and var1(8-1,*) eq step and var1(10-1,*) eq nfor and var1(4-1,*) eq nk)
    ind2 = where(var2(0,*) eq var3zl(nn-1) and var2(2-1,*) eq levtype and var2(8-1,*) eq step and var2(10-1,*) eq nfor and var2(4-1,*) eq nk)
   ;print, var3zl(nn-1), ' L', nk, ' exp',expref, ' bias:', var2(16-1,ind2), 'rms:', var2(18-1,ind2)
    bias = abs(float(var1(16-1,ind1))) / ( abs(float(var1(16-1,ind1))) + abs(float(var2(16-1,ind2))) )
    rms  =     float(var1(18-1,ind1))  / (     float(var1(18-1,ind1))  +     float(var2(18-1,ind2))  )
    if ( abs(float(var1(16-1,ind1))) + abs(float(var2(16-1,ind2))) ) lt 1e-8 then begin
      bias = 0.5
    endif
    if (     float(var1(18-1,ind1))  +     float(var2(18-1,ind2))  ) lt 1e-8 then begin
      rms  = 0.5
    endif
  
    xxx =-0.11+nn*0.14
    yyy =0.93
    vert=0.01
  
    xyouts, xxx+0.03, yyy+0.015, var3zl(nn-1), /normal, charsize=0.7
    header = 'level    rms      new'
    xyouts, xxx, yyy, header, /normal, charsize=0.7
    err1=var1(18-1,ind1)
    err2=var2(18-1,ind2)
    err =rms
    xoff=0.0
    amp =40.0
    scale = (err-0.5) * (err gt 0.5)
    xr= [0.0,0.01,0.01,0.0 ,0.0] * amp * scale(0) + xxx + xoff
    yr= [0.0,0.0 ,0.01,0.01,0.0] * 0.8            + yyy
    polyfill, xr   +0.04, yr   -nk*vert,              /normal, color=240
    scale=(0.5-err) * (err lt 0.5)
    xr=-[0.0,0.01,0.01,0.0 ,0.0] * amp * scale(0) + xxx + xoff
    polyfill, xr     +0.04, yr     -nk*vert,           /normal, color=160
    plots,    xr(3:4)+0.04, [yr(3)+0.002,yr(4)-0.002]-nk*vert,           /normal

    xyouts,   xr(0)     , yr(0)-nk*vert, strcompress(string(nk)), /normal, charsize=0.5
    xyouts,   xr(0)+0.07, yr(0)-nk*vert, err1,  /normal, charsize=0.4
   end
  end
  
  erase
  
;--------------------------------------------------------------------------------------------

  levtype='sfc'

  xyouts, 0.15, 0.97, /normal, charsize=1.0, $
    'Comparison of exp'+expnum+' to exp'+expref+ $
    '   forecast hours: '+step+'   number of forecasts: '+nfor+ ' from '+inidate+'   '+levtype
  
  for nn=1,nvar2d do begin
   ;print,''
   ;print,'searching for variable ', var2d(nn-1)
    ind1 = where(var1(0,*) eq var2d(nn-1) and var1(2-1,*) eq levtype and var1(8-1,*) eq step and var1(10-1,*) eq nfor)
    ind2 = where(var2(0,*) eq var2d(nn-1) and var2(2-1,*) eq levtype and var2(8-1,*) eq step and var2(10-1,*) eq nfor)
   ;print, var2d(nn-1), ' exp',expref, ' bias:', var2(16-1,ind2), 'rms:', var2(18-1,ind2)  
    bias = abs(float(var1(16-1,ind1))) / ( abs(float(var1(16-1,ind1))) + abs(float(var2(16-1,ind2))) )
    rms  =     float(var1(18-1,ind1))  / (     float(var1(18-1,ind1))  +     float(var2(18-1,ind2))  )
  
    xxx=0.1
    yyy=0.9
    vert=0.033
  
    header = 'variable                    bias             new,ref            rms             new,ref'
    xyouts, xxx, yyy, header, /normal, charsize=0.8
    for nerr=1,2 do begin
      if nerr eq 1 then begin
        err1=var1(16-1,ind1)
        err2=var2(16-1,ind2)
        err =bias
        xoff=0.0
        amp =12.0
      endif else begin
        err1=var1(18-1,ind1)
        err2=var2(18-1,ind2)
        err =rms
        xoff=0.14
        amp =30.0
      endelse
      scale = (err-0.5) * (err gt 0.5)
      xr= [0.0,0.01,0.01,0.0 ,0.0] * amp * scale(0) + xxx + xoff
      yr= [0.0,0.0 ,0.01,0.01,0.0] * 1.0            + yyy
      polyfill, xr   +0.12 , yr   -nn*vert,              /normal, color=240
      scale=(0.5-err) * (err lt 0.5)
      xr=-[0.0,0.01,0.01,0.0 ,0.0] * amp * scale(0) + xxx + xoff
      polyfill, xr     +0.12, yr  -nn*vert,           /normal, color=160
      plots,    xr(3:4)+0.12, [yr(3)+0.002,yr(4)-0.002]-nn*vert,           /normal

      if nerr eq 1 then begin
      xyouts,   xr(0)     , yr(0)-nn*vert, var2d(nn-1), /normal, charsize=0.7
      endif
      xyouts,   xr(0)+0.18, yr(0)-nn*vert+0.006, err1,  /normal, charsize=0.5
      xyouts,   xr(0)+0.18, yr(0)-nn*vert-0.006, err2,  /normal, charsize=0.5
    end
   
  end

  erase

end

paperclose


end
