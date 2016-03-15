; batch execution of ICON scores plotting
;-----------------------------------------------------------------------
; Useage: run as: idl batch_scores -args [direxp] [dirref] [expnum]
; [expref] [inidate] [step] [nfor]
;
; e.g.:
; idl batch_scores -args /fe1-daten/mkoehler/plots/icon 82 81 20110101
;
; Martin Koehler, Sep 2012 & Dec 2012
;-----------------------------------------------------------------------

 args  = COMMAND_LINE_ARGS(count=nargs)
 direxp  = args(0)          ;long(args(0)), fix(args(2))
 dirref  = args(1)
 expnum  = args(2)
 expref  = args(3)
 inidate = args(4)
 step    = args(5)
 nfor    = args(6)

 print,'Input arguments: direxp=    ',direxp
 print,'                 dirref=    ',dirref
 print,'                 expnum=    ',expnum
 print,'                 expref=    ',expref
 print,'                 inidate=   ',inidate
 print,'                 step=      ',step
 print,'                 nfor=      ',nfor

 .run scores

 outdate   = long(inidate) + step/24
 plotfile2 = strcompress(string(outdate),/remove_all)+'00.L1.ps'
 
; mean plots

 plotfile1 = direxp+'/scores.NWP.r2B06L90.'+expnum+'_vs_'+expref

 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'ml', 'mean'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.ml.'  + plotfile2
 
 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'pl', 'mean'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.pl.'  + plotfile2

 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'zl', 'mean'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.zl.'  + plotfile2

 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'sfc', 'mean'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.sfc.' + plotfile2

; RMS plots

 plotfile1 = direxp+'/scores.NWP.r2B06L90.'+expnum+'_rms_vs_'+expref
 
 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'ml', 'rms'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.ml.'  + plotfile2
 
 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'pl', 'rms'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.pl.'  + plotfile2

 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'zl', 'rms'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.zl.'  + plotfile2

 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'sfc', 'rms'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.sfc.' + plotfile2


 exit

