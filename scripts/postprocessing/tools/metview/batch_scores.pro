; batch execution of ICON scores plotting
;-----------------------------------------------------------------------
; Useage: run as: idl batch_scores -args  [plotdir] [expnum] [expref] [inidate]
; e.g.:
; idl batch_scores -args /fe1-daten/mkoehler/plots/icon 82 81 20110101
;
; Martin Koehler, Sep 2012
;-----------------------------------------------------------------------

 args  = COMMAND_LINE_ARGS(count=nargs)
 plotdir = args(0)          ;long(args(0)), fix(args(2))
 expnum  = args(1)
 expref  = args(2)
 inidate = args(3)

 print,'Input arguments: plotdir=   ',plotdir
 print,'                 expnum=    ',expnum
 print,'                 expref=    ',expref
 print,'                 inidate=   ',inidate

.run scores
 scores, plotdir, expnum, expref, inidate

 spawn,'\cp -f plot.ps '+plotdir+"/nwp.exp"+expnum+"/scores_exp"+expnum+".ps"

 exit
