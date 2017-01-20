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
        parser.add_argument("startdate",      help="start date [YYYYMMDDhhmmss]")
        parser.add_argument("enddate",        help="end date [YYYYMMDDhhmmss]")
        parser.add_argument("--increment",    help="increment time span [h]",                                 default=24)
        parser.add_argument("--mode",         help="mode: "+str(MODE_IAU)+"=IAU, "+str(MODE_NOIAU)+"=No IAU", default=MODE_IAU)
        parser.add_argument("--ensemble",     help="read ensemble data",  dest="ensemble", action='store_true')
        parser.add_argument("--emember",      help="ensemble member",                                         default="1")
        parser.set_defaults(ensemble=False)
        args = parser.parse_args()

        cur_date  = datetime.datetime.strptime(args.startdate, DATEFMT)
        enddate   = datetime.datetime.strptime(args.enddate,   DATEFMT)
        increment = datetime.timedelta(hours=int(args.increment))
        delta_3h  = datetime.timedelta(hours=3)

        # ------------------------------------------------------------
        # loop over time samples:

        filelist = ""
        while (cur_date <= enddate):
            # format current date as string
            datestr      = datetime.datetime.strftime(cur_date, DATEFMT)
            datestr_3    = datetime.datetime.strftime(cur_date - delta_3h, DATEFMT) # = "cur_date - 3h"
            print "date: ", datestr
            # calculate current day, 00UTC:
            cur_date_00 = cur_date.replace(hour=0, minute=0, second=0, microsecond=0)

            # ------------------------------------------------------------
            # depending on ensemble mode and day hour: make some
            # special settings
            if (not args.ensemble):
                # non-ensemble mode
                cat       = "icogl130l90_main_an_rout"
                cat_eu    = "icoeu_main_an_rout" 
                cat_fg    = "icogl130l90_pre_fc_rout"
                cat_fg_eu = "icoeu_pre_fc_rout" 
                # lin=20 (localInformationNumber) is related to the
                # W_SO analysis including SMA.
                lin       = "lin=20"
            elif (args.ensemble):
                # ensemble mode
                cat       = "icogle_main_an_rout letype=101 enum="+args.emember
                cat_eu    = "icoeue_main_an_rout letype=101 enum="+args.emember
                cat_fg    = "icogle_pre_fc_rout  letype=101 enum="+args.emember
                cat_fg_eu = "icoeue_pre_fc_rout  letype=101 enum="+args.emember
                # lin=0 (localInformationNumber): In case of W_SO
                # variable contains the ensemble perturbations (?) but
                # no SMA increments.
                lin       = "lin=0"
            else:
                raise Exception("Unknown resolution!")

            read_w_so = True
            if (cur_date == cur_date_00):
                # if hour == 0 then gptype_str=0, else -1 for
                # t_so,t_sea,fr_ice,h_ice,t_ice:
                gptype_str = "0"
            else:
                gptype_str = "-1"
                # if hour /= 0 and not ensemble then do not read w_so
                if (not ensemble):
                    read_w_so = False

            # ------------------------------------------------------------
            # assemble "pbank" requests:
            #  gptype  ===  typeOfGeneratingProcess
            #  lin     ===  localInformationNumber
            #
            req_filename = os.environ['TMPDIR']+"/skyreq_"+datestr+"00R_G+N"
            f = open(req_filename, "w")
            f.write("reqColl proc=parallel timeout=0\n")
            prefix    = "read db=roma bin info=countPlus "

            if (args.mode == MODE_IAU):
                f.write((prefix + "cat="+cat+"       d={0} gptype=201   f=igaf{0}R.grb  p=u,v,t,p,qv                                \n").format(datestr))
                f.write((prefix + "cat="+cat+"       d={0} gptype={1}   f=igaf{0}R.grb  p=t_so,t_sea,fr_ice,h_ice,t_ice             \n").format(datestr,gptype_str))
                if (read_w_so):
                    f.write((prefix + "cat="+cat+"   d={0} gptype=201   f=igaf{0}R.grb  p=w_so "+lin+"                              \n").format(datestr))
                f.write((prefix + "cat="+cat+"       d={0} gptype=201   f=igaf{0}R.grb  p=freshsnw,h_snow                           \n").format(datestr))
                f.write((prefix + "cat="+cat_fg+"    d={0} s[s]=5400    f=igfff{0}-0130R.grb  lvt1=!100                             \n").format(datestr_3))
                f.write((prefix + "cat="+cat_eu+"    d={0} gptype=0     f=ieaf{0}R.grb  p=t_so,t_sea,fr_ice,h_ice,t_ice             \n").format(datestr))
                f.write((prefix + "cat="+cat_eu+"    d={0} gptype=201   f=ieaf{0}R.grb  p=w_so "+lin+"                              \n").format(datestr))
                f.write((prefix + "cat="+cat_eu+"    d={0} gptype=201   f=ieaf{0}R.grb  p=freshsnw,h_snow                           \n").format(datestr))
                f.write((prefix + "cat="+cat_fg_eu+" d={0} s[s]=5400    f=iefff{0}-0130R.grb  lvt1=!100                             \n").format(datestr_3))
            elif (args.mode == MODE_NOIAU):
                f.write((prefix + "cat="+cat+"       d={0} gptype=!201  f=igaf{0}R.grb p=u,v,t,p,qv,hhl                             \n").format(datestr))
                f.write((prefix + "cat="+cat+"       d={0} gptype={1}   f=igaf{0}R.grb p=t_so,t_sea,fr_ice,h_ice,t_ice              \n").format(datestr, gptype_str))
                if (read_w_so):
                    f.write((prefix + "cat="+cat+"   d={0} gptype=0     f=igaf{0}R.grb p=w_so "+lin+"                               \n").format(datestr))
                f.write((prefix + "cat="+cat+"       d={0} gptype=0     f=igaf{0}R.grb p=w_snow,t_snow,w_i,rho_snow,freshsnw,h_snow \n").format(datestr))
                f.write((prefix + "cat="+cat_fg+"    d={0} s[s]=10800   f=igfff{0}-0300R.grb lvt1=!100                              \n").format(datestr_3))
                f.write((prefix + "cat="+cat_eu+"    d={0} gptype=0     f=ieaf{0}R.grb p=t_so,t_sea,fr_ice,h_ice,t_ice,hhl          \n").format(datestr))
                f.write((prefix + "cat="+cat_eu+"    d={0} gptype=0     f=ieaf{0}R.grb p=w_so "+lin+"                               \n").format(datestr))
                f.write((prefix + "cat="+cat_eu+"    d={0} gptype=0     f=ieaf{0}R.grb p=w_snow,t_snow,w_i,rho_snow,freshsnw,h_snow \n").format(datestr))
                f.write((prefix + "cat="+cat_fg_eu+" d={0} s[s]=10800   f=iefff{0}-0300R.grb lvt1=!100                              \n").format(datestr_3))
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
