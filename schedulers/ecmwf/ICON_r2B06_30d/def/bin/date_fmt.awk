# DATE_FMT.AWK
#
# reformatting of a date string
# from format "2011010100"
# to format "2011-01-01T00:00:00Z"
#
# Usage:
#  awk -v in_date="2011010100" -f ../bin/date_fmt.awk 
#
# Corresponding author:
#   Florian Prill, DWD, mailto:florian.prill@dwd.de
BEGIN {
  YY = substr(in_date,1,4);
  MM = substr(in_date,5,2);
  DD = substr(in_date,7,2);
  HH = substr(in_date,9,2);
  MM = 00
  SS = 00
  printf("%04d-%02d-%02dT%02d:%02d:%02dZ", YY, MM, DD, HH, MM, SS);
}
