#!/usr/bin/env python

# SKY4ICON.PY
# Retrieve ICON data from the DWD "Sky" database.
#
# 01/2017 : D. Reinert/F. Prill, DWD

import argparse, datetime, os, subprocess


# ------------------------------------------------------------
# CONSTANTS

DATEFMT    = '%Y%m%d%H%M%S'
MODE_IAU   = 1
MODE_NOIAU = 2

# ------------------------------------------------------------
# Main program (wrapped by a subroutine)
#
def main():

    print """\nSKY4ICON.PY
    Retrieve ICON data from the DWD "Sky" database.
    01/2017 : D. Reinert/F. Prill, DWD\n"""
    
    try:
    
        # parse command-line options
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("startdate",   help="start date [YYYYMMDDhhmmss]")
        parser.add_argument("enddate",     help="end date [YYYYMMDDhhmmss]")
        parser.add_argument("--increment", help="increment time span [h]", default=24)
        parser.add_argument("--mode",      help="mode: "+str(MODE_IAU)+"=IAU, "+str(MODE_NOIAU)+"=No IAU", default=MODE_IAU)
        args = parser.parse_args()

        cur_date  = datetime.datetime.strptime(args.startdate, DATEFMT)
        enddate   = datetime.datetime.strptime(args.enddate,   DATEFMT)
        increment = datetime.timedelta(hours=int(args.increment))
        delta_3h  = datetime.timedelta(hours=3)

        # loop over time samples:
        filelist = ""
        while (cur_date < enddate):
            datestr      = datetime.datetime.strftime(cur_date, DATEFMT)
            datestr_3    = datetime.datetime.strftime(cur_date - delta_3h, DATEFMT) # = "cur_date - 3h"

            # create a "pbank" request
            req_filename = os.environ['TMPDIR']+"/skyreq_"+datestr+"00R_G+N"
            f = open(req_filename, "w")
            f.write("reqColl proc=parallel timeout=0")
            prefix = "read db=roma bin info=countPlus "
            if (args.mode == MODE_IAU):
                f.write(prefix + "cat=icogl130l90_main_an_rout d={0} gptype=201   f=igaf{0}R.grb  p=u,v,t,p,qv                \n".format(datestr))
                f.write(prefix + "cat=icogl130l90_main_an_rout d={0} gptype=0     f=igaf{0}R.grb  p=t_so,fr_ice,h_ice,t_ice   \n".format(datestr))
                f.write(prefix + "cat=icogl130l90_main_an_rout d={0} gptype=201   f=igaf{0}R.grb  p=w_so lin=20               \n".format(datestr))
                f.write(prefix + "cat=icogl130l90_main_an_rout d={0} gptype=201   f=igaf{0}R.grb  p=freshsnw,h_snow           \n".format(datestr))
                f.write(prefix + "cat=icogl130l90_pre_fc_rout  d={0} s[s]=5400    f=igfff{0}-0130R.grb  lvt1=!100             \n".format(datestr_3))
                f.write(prefix + "cat=icoeu_main_an_rout       d={0} gptype=0     f=ieaf{0}R.grb  p=t_so,fr_ice,h_ice,t_ice   \n".format(datestr))
                f.write(prefix + "cat=icoeu_main_an_rout       d={0} gptype=201   f=ieaf{0}R.grb  p=w_so lin=20               \n".format(datestr))
                f.write(prefix + "cat=icoeu_main_an_rout       d={0} gptype=201   f=ieaf{0}R.grb  p=freshsnw,h_snow           \n".format(datestr))
                f.write(prefix + "cat=icoeu_pre_fc_rout        d={0} s[s]=5400    f=iefff{0}-0130R.grb  lvt1=!100             \n".format(datestr_3))
            elif (args.mode == MODE_NOIAU):
                f.write(prefix + "cat=icogl130l90_main_an_rout d={0} gptype=!201  f=igaf{0}R.grb p=u,v,t,p,qv,hhl             \n".format(datestr))
                f.write(prefix + "cat=icogl130l90_main_an_rout d={0} gptype=0     f=igaf{0}R.grb p=t_so,fr_ice,h_ice,t_ice    \n".format(datestr))
                f.write(prefix + "cat=icogl130l90_main_an_rout d={0} gptype=0     f=igaf{0}R.grb p=w_so lin=20                \n".format(datestr))
                f.write(prefix + "cat=icogl130l90_main_an_rout d={0} gptype=0     f=igaf{0}R.grb p=w_snow,t_snow,w_i,rho_snow,freshsnw,h_snow\n".format(datestr))
                f.write(prefix + "cat=icogl130l90_pre_fc_rout  d={0} s[s]=10800   f=igfff{0}-0300R.grb lvt1=!100              \n".format(datestr_3))
                f.write(prefix + "cat=icoeu_main_an_rout       d={0} gptype=0     f=ieaf{0}R.grb p=t_so,fr_ice,h_ice,t_ice,hhl\n".format(datestr))
                f.write(prefix + "cat=icoeu_main_an_rout       d={0} gptype=0     f=ieaf{0}R.grb p=w_so lin=20                \n".format(datestr))
                f.write(prefix + "cat=icoeu_main_an_rout       d={0} gptype=0     f=ieaf{0}R.grb p=w_snow,t_snow,w_i,rho_snow,freshsnw,h_snow\n".format(datestr))
                f.write(prefix + "cat=icoeu_pre_fc_rout        d={0} s[s]=10800   f=iefff{0}-0300R.grb lvt1=!100              \n".format(datestr_3))
            else:
                raise Exception("Unknown mode!")

            f.close()
            filelist += req_filename + " ";
            cur_date += increment

        subprocess.call(["pbank", filelist])  # execute "pbank"
        print "    * " + datetime.datetime.now().strftime("%b %d %Y %H:%M:%S") + " :: done."    
    
    except Exception as e:
        print str(e)



# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------

if __name__ == '__main__':
    main()
