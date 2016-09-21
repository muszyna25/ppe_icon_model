# BUILD_VARLISTS
#
# Python script generating variable lists from LaTeX tables.
#
# Initial revision:     2015-04-16 : F. Prill, DWD
# Add lists for int2lm: 2015-04-22 : H. Frank, DWD

import re, os, commands, itertools
import argparse

# parse command line arguments (LaTeX filename)
parser = argparse.ArgumentParser()
parser.add_argument("file", type=str, help="input LaTeX file name")
args = parser.parse_args()

# result lists:
list_tri_glb = []
list_tri_loc = []
list_tri_i2l = []
list_ll_glb  = []
list_ll_loc  = []

re_groups  = re.compile(r'^\s*\\groups')
re_column1 = re.compile(r'(^\s*\\onlyglb\s*{)?(^\s*\\onlyloc\s*{)?(\s*\\groups\[)([^\]]*)(\]\[)([^\]]*)(\])')
re_tri_ll  = re.compile(r'(^\s*\\onlyglb\s*{)?(^\s*\\onlyloc\s*{)?(\s*tri|ll\s*)(\})?')
re_i2l     = re.compile(r'\si2l\s*$')
re_pl      = re.compile(r'\spl\s*$')
re_comment = re.compile(r'^\s*%')

print "> Open '" + args.file.strip() + "' and scan for output field short names..."
with open(args.file) as f:
    for line in f.readlines():
        reg = re.search(re_comment, line)
        if ( reg != None ):
            continue
        reg = re.search(re_column1, line)
        columns = line.split('&')
#       if len(columns)>3:
        if ( reg != None and len(columns)>3):

            # ignore comment lines

            # short name: replace LaTeX "\_" and remove trailing LaTeX commands
            # like "\footnotemark..."
            shortname = columns[1].replace('\_','_').split('\\')[0].strip()

            # test pl marker. Exclude data on pressure levels in the lists
            regexpl = re.search(re_pl, columns[-1])
            if ( regexpl != None):
                continue

            # regular expression for the first table column
#           re_column1 = re.compile(r'(^\s*\\onlyglb\s*{)?(^\s*\\onlyloc\s*{)?(\s*\\groups\[)([^\]]*)(\]\[)([^\]]*)(\])')
            regexp = re.search(re_column1, columns[0])
            all_glb = False
            if (regexp.group(1) != None):
                if (regexp.group(1) != None):
                    all_glb = True
            all_loc = False
            if (regexp.group(2) != None):
                if (regexp.group(2) != None):
                    all_loc = True
            # test "tri" marker:
#           re_tri_ll = re.compile(r'(^\s*\\onlyglb\s*{)?(^\s*\\onlyloc\s*{)?(\s*tri|ll\s*)(\})?')
            regexp2 = re.search(re_tri_ll, str(regexp.group(4)))
            tri_loc = False
            tri_glb = False
            if (regexp2 != None):
                tri_glb = (regexp2.group(2) == None)
                tri_loc = (regexp2.group(1) == None)
            tri_loc &= not all_glb
            tri_glb &= not all_loc
            # test "ll" marker:
            regexp2 = re.search(re_tri_ll, str(regexp.group(6)))
            ll_loc = False
            ll_glb = False
            if (regexp2 != None):
                ll_glb = (regexp2.group(2) == None)
                ll_loc = (regexp2.group(1) == None)
            ll_loc &= not all_glb
            ll_glb &= not all_loc

            # test i2l marker
            regexp3 = re.search(re_i2l, columns[-1])
            i2l = False
            if ( regexp3 != None):
                i2l = True
	    tri_i2l = False
            if ( i2l):
		    if ( tri_glb):
			    tri_i2l = True
			    tri_glb = False
			    tri_loc = False

            if (tri_i2l == True):  list_tri_i2l.append(shortname)
            if (tri_glb == True):  list_tri_glb.append(shortname)
            if (tri_loc == True):  list_tri_loc.append(shortname)
            if (ll_glb == True):   list_ll_glb.append(shortname)
            if (ll_loc == True):   list_ll_loc.append(shortname)

            #print columns[0]
            #print shortname, ": tri glb/loc = ", tri_glb, tri_loc, "; ll glb/loc = ", ll_glb, ll_loc

print "\n> list of variables for int2lm, triangular grid:"
print list_tri_i2l
print "\n> list of variables for triangular (not int2lm), global grid):"
print list_tri_glb
print "\n> list of variables for triangular (not int2lm), local grid:"
print list_tri_loc
print "\n> list of variables for lon-lat, global grid:"
print list_ll_glb
print "\n> list of variables for lon-lat, local grid:"
print list_ll_loc
