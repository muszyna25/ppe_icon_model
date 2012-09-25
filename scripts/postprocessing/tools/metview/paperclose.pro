 pro paperclose

 device, /close_file, /helvetica
;spawn,'lpr plot.ps'
;spawn,'rm -f plot.ps'
 set_plot,'x'

 !p.charsize=1 & !p.charthick=1 & !p.font=-1
 !p.thick=1 & !x.thick=1 & !y.thick=1
 !p.position=0
;!x.margin = 0 & !y.margin = 0
 !p.region = 0

 end
