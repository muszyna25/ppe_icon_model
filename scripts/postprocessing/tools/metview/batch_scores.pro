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

 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'ml'
 spawn,'\cp -f plot.ps metplots/scores_'+expnum+'_vs_'+expref+'_'+ $
   inidate+'+'+step+'_ml.ps'

 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'pl'
 spawn,'\cp -f plot.ps metplots/scores_'+expnum+'_vs_'+expref+'_'+ $
   inidate+'+'+step+'_pl.ps'

 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'zl'
 spawn,'\cp -f plot.ps metplots/scores_'+expnum+'_vs_'+expref+'_'+ $
   inidate+'+'+step+'_zl.ps'

 scores, direxp, dirref, expnum, expref, inidate, step, nfor, 'sfc'
 spawn,'\cp -f plot.ps metplots/scores_'+expnum+'_vs_'+expref+'_'+ $
   inidate+'+'+step+'_sfc.ps'

 exit
