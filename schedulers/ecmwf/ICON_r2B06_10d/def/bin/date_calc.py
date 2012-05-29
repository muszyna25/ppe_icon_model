import optparse
import datetime

#  -------------------------------------------------
#  REQUEST_DATE.PY
#  -------------------------------------------------
#
#  Python script doing some date calculations.
#
#  This script is part of the ICON SMS suite
#  Initial implementation: F. Prill, DWD (2012-05-07)
#
#  Corresponding author:
#    Florian Prill, DWD, mailto:florian.prill@dwd.de

def main():
    p = optparse.OptionParser()
    p.add_option('--date',     '-d', default="2012050100",  help="requested date YYYMMDDHH")
    p.add_option('--subtract', '-s', default="0",           help="days to subtract")
    p.add_option('--action',   '-a', default="printdate",   help="requested type of output: printdate/weekday/printfmt")
    #
    options, arguments = p.parse_args()

    # parse time string:
    date_object  = datetime.datetime.strptime(options.date, '%Y%m%d%H')
    # date calculation:
    date_object2 = date_object - datetime.timedelta(int(options.subtract))

    if (options.action == "printdate"):
        # print date (formatted output)
        print(date_object2.strftime("%Y%m%d%H")) 
    elif (options.action == "weekday"):
        # print only day of week
        print(date_object2.strftime("%w")) 
    elif (options.action == "printfmt"):
        # format as required by ICON:
        print(date_object2.strftime("%Y-%m-%dT%H:00:00Z"))
    else:
        print "Unknown command line option \"", options.action,"\"!"


if __name__ == '__main__':
    main()

