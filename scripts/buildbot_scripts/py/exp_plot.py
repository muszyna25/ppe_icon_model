#================================== NEW ================================

from buildbot.status.web.base import HtmlResource
from twisted.web.util import Redirect
import os, time
from buildbot import version, util
from datetime import datetime, timedelta

class EXP_plot(HtmlResource):
    title = "Exp Plots"

    def head(self, request):
        h='''
<link href="archive.css" rel="stylesheet" type="text/css" />

<script type='text/javascript'>
 function anzeigen(das,was){
 if (document.getElementById(das).style.display=='none') {
  document.getElementById(das).style.display='block'; 
  was.src="open.gif";
 }
 else {
  document.getElementById(das).style.display='none';
   was.src='closed.gif';
  }
}
 </script>
        '''
        return h
						
   
    def body(self, request):
        global exp_plot_Info
	global exp_file_Info
        global l_nightly
	self.get_info(request)
        status = self.getStatus(request)	
        Archive_Button_Dict = {}

# -----------------------------------------------------------
        def add_file(p,e,Dict):
	  ret = False
	  file_L = []
	  p += e + "/plots/"
	  for f in os.listdir(p):
	    file_L.append(f)
       	    Archive_Button_Dict['file'].append(f)
	    ret = True
	    
# If ret == True we got a file in the list
	  if ret:
            Dict[e] = file_L
	    
          return ret

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
	def add_exp(p,b,Dict):
	  ret = False
	  exp_D = {}
	  p += b + "/"
	  for e in os.listdir(p):
	    if e.find("test_") == 0:
	      if e == exp_plot_Info or exp_plot_Info == "all":
	        if add_file(p,e,exp_D):
		  Archive_Button_Dict['exp'].append(e)
	          ret = True

# If ret == True we found an existing experiment mit a file list
	  if ret:
            Dict[b] = exp_D
          
	  return ret

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
	def add_build(p,c,Dict):
	  ret = False
	  build_D = {}
	  p += c + "/"
	  for b in os.listdir(p):
	    if add_exp(p,b,build_D):
	      Archive_Button_Dict['build'].append(b)
	      ret = True

# If ret == True we found a build with correct sub-information
	  if ret:
            Dict[c] = build_D
          
	  return ret

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
	def add_comp(p,r,r_Dict):
	  ret = False
	  comp_D = {}
	  p += r + "/"
	  for c in os.listdir(p):
	    if add_build(p,c,comp_D):
  	      Archive_Button_Dict['comp'].append(c)
	      ret = True

# If ret == True we found a computer with correct sub-information
	  if ret:
            r_Dict[r] = comp_D
          
          return ret

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
	def add_rev(p,D,d_Dict):
	  ret = False
	  rev_D = {}
	  p += D + "/buildbot/"
	  for r in os.listdir(p):
            if (r >= svn_Rev_from) and (r <= svn_Rev_to):
	      if add_comp(p,r,rev_D):
	        Archive_Button_Dict['rev'].append(r)
	        ret = True
#	    else:
#	      print "Rev: "+r+" not used"

# If ret == True we found a revision number with correct sub-information
	  if ret:
            d_Dict[D] = rev_D
          
          return ret
	      
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        def create_Exp_Dict(d,r,c,b,eD):
	  build_D = {}
	  comp_D  = {}
	  rev_D   = {}
	  data_D  = {}
	  build_D[b] = eD
	  comp_D[c] = build_D
	  rev_D[r] = comp_D
	  data_D[d] = rev_D
	  return data_D 
	  
        def create_data_Dict(r,c,b,eD):
	  build_D = {}
	  comp_D  = {}
	  rev_D   = {}
	  build_D[b] = eD
	  comp_D[c] = build_D
	  rev_D[r] = comp_D
	  return rev_D 
        
	def create_rev_Dict(c,b,eD):
	  build_D = {}
	  comp_D  = {}
	  build_D[b] = eD
	  comp_D[c] = build_D
	  return comp_D 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	def create_comp_Dict(b,eD):
	  build_D = {}
	  build_D[b] = eD
	  return build_D 
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	def Exp_add(Exp_Dict,d,r,c,b,eD):
	   
	   for e in eD:
             if e in Exp_Dict:
	       if d in Exp_Dict[e]:
	         print "date existiert " + d
	         if r in Exp_Dict[e][d]:
	           print "ref existiert " + r
	           if c in Exp_Dict[e][d][r]:
	             print "comp existiert " + c
	             if b in Exp_Dict[e][d][r][c]:
	               print "build existiert " + c
	               Exp_Dict[e][d][r][c][b] = eD[e]
		     else:
	               Exp_Dict[e][d][r][c][b] = eD[e]
	           else:
	             Exp_Dict[e][d][r][c] = create_comp_Dict(b,eD[e])
	         else:
	           Exp_Dict[e][d][r] = create_rev_Dict(c,b,eD[e])
	       else:
	         Exp_Dict[e][d] = create_data_Dict(r,c,b,eD[e])
             else:
	       Exp_Dict[e] = create_Exp_Dict(d,r,c,b,eD[e])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
        def Comp_create_date(d,r,eD):
	  rev_D   = {}
	  data_D  = {}
	  rev_D[r] = eD
	  data_D[d] = rev_D
	  return data_D 
	
        def Comp_create_rev(r,eD):
	  rev_D   = {}
	  rev_D[r] = eD
	  return rev_D 
	
	def Comp_add(Comp_Dict,d,r,cD):
	  for c in cD:
           if c in Comp_Dict:
             print "comp exists " + c
             if d in Comp_Dict[c]:
               print "date exists " + d
               if r in Comp_Dict[c][d]:
                 print "!!!! rev exists " + r
                 Comp_Dict[c][d][r] = cD[c]
	       else:
                 Comp_Dict[c][d][r] = cD[c]
	     else:
	       Comp_Dict[c][d] = Comp_create_rev(r,cD[c])
	   else:
	     Comp_Dict[c] = Comp_create_date(d,r,cD[c])
# -----------------------------------------------------------
	
	if exp_plot_Info == "NotSet":
	  data = self.footer(status,request)
          return data

	Archive_Button_Dict['date']  = []
	Archive_Button_Dict['rev']   = []
	Archive_Button_Dict['comp']  = []
	Archive_Button_Dict['build'] = []
	Archive_Button_Dict['exp']   = []
	Archive_Button_Dict['file']  = []
	
#        ExpListInfo = []

	
        if l_nightly:
          file = open("icon_run/svn_info")
          while 1:
            line = file.readline()
            if not line:
              break
          
            if line.find("Revision:") >= 0:
	      tmp,tmp_Rev = line.split("Revision:",1)
              svn_Rev_from = tmp_Rev.strip(" ").replace("\n","")
	      svn_Rev_to = svn_Rev_from
 
	      svn_Date_to = time.strftime("%Y-%m-%d",time.localtime(util.now()))
	      now_time = time.strftime("%H%M",time.localtime(util.now()))
	      
	      yesterday = datetime.now() - timedelta(days=1)
	      
	      if now_time > "2200":
		svn_Date_from = svn_Date_to
	      else:
	        svn_Date_from = yesterday.strftime("%Y-%m-%d")
		      
	      
              print "====== svn_Date_from_to ======"
	      print now_time
	      print svn_Date_from
	      print svn_Date_to
	      
              print "====== svn_Date_from_to ======"
	      
#            if line.find("Last Changed Date:") >= 0:
#	      tmp,tmp_Date = line.split("Last Changed Date:",1)
#              tmp_Date.replace("\n","")
#              tmp_Date.strip(" ")
#              svn_Date = tmp_Date.split(' ')[1].replace("\n","")
               
          file.close
	else:
          svn_Rev_from  = "9440"
          svn_Rev_to    = "9450"
#	  svn_Date_from = "2012-06-01"
#	  svn_Date_to   = "2012-06-31"

#========================================================================================================
#
# Begin searching for plots
#
#========================================================================================================

# Build new emty date_Dict type

        

# Build new emty rev_Dict type and save it in the date_Dict under the correct Date
	date_Dict = {}
        
	p_date = "public_html/archive/"
        print "==== WS ===="
	for DATE in os.listdir(p_date):
	  if DATE.find("2012-06") == 0:
            if (DATE >= svn_Date_from) and (DATE <= svn_Date_to):
	      if add_rev(p_date,DATE,date_Dict):
	        Archive_Button_Dict['date'].append(DATE)
#	      Archive_Button_Dict['date'].append(DATE)
#	      print "    Used:" + DATE
#	    else:
#	      print "Not Used:" + DATE
	    
# Build a new List where the Experiment Name ist the hightest index
	Exp_Dict = {}
	for d in date_Dict:
	  for r in date_Dict[d]:
	    for c in date_Dict[d][r]:
	      for b in date_Dict[d][r][c]:
	        Exp_add(Exp_Dict,d,r,c,b,date_Dict[d][r][c][b])
             
# Build a new List where the Experiment Name ist the hightest index
	Comp_Dict = {}
	for d in date_Dict:
	  for r in date_Dict[d]:
	     Comp_add(Comp_Dict,d,r,date_Dict[d][r])

	print "====== Exp_add ======"
        print "Exp_Dict:"
        print Exp_Dict
	print "====== Exp_add ======"
 
	print "====== Comp_add ======"
        print "Comp_Dict:"
        print Comp_Dict
	print "====== Exp_add ======"
        Archive_Button_Dict['date']  = list(sorted(set(Archive_Button_Dict['date'])))
 	Archive_Button_Dict['rev']   = list(sorted(set(Archive_Button_Dict['rev'])))
	Archive_Button_Dict['comp']  = list(sorted(set(Archive_Button_Dict['comp'])))
	Archive_Button_Dict['build'] = list(sorted(set(Archive_Button_Dict['build'])))
	Archive_Button_Dict['exp']   = list(sorted(set(Archive_Button_Dict['exp'])))
	Archive_Button_Dict['file']  = list(sorted(set(Archive_Button_Dict['file'])))
        
#        print "Archive_Button_Dict"
#	print Archive_Button_Dict['date']
# 	print Archive_Button_Dict['rev']
#	print Archive_Button_Dict['comp']
#	print Archive_Button_Dict['build']
#	print Archive_Button_Dict['exp']
#       print Archive_Button_Dict['file']
#        print "date_Dict"
#        print date_Dict
#	print "==== WS ===="
	
#	if len(Archive_Button_Dict['comp']) > 1:
#	  Archive_Button_Dict['comp'].insert(0, 'all')
	  
#	if len(Archive_Button_Dict['build']) > 1:
#	  Archive_Button_Dict['build'].insert(0, 'all')
#	if len(Archive_Button_Dict['exp']) > 1:
#	  Archive_Button_Dict['exp'].insert(0, 'all')
#	if len(Archive_Button_Dict['file']) > 1:
#	  Archive_Button_Dict['file'].insert(0, 'all')

#========================================================================================================
#
# End searching for plots
#
#========================================================================================================

#==============================================
# Begin of whole div
        
	data = "<div id=\"page_menu\" >\n"

# Begin of left side  div
	
	data += "<div style=\"float:left; padding:3px; margin:5px width:350px;\">\n"

# Begin of ´ICON Buildbot Archive´ div
        
        data += "<div id=\"arch_menu\" >\n"
	data += "<h1>ICON Buildbot Archive</h1>\n"
#        data += "<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
#	data += "<tr><td style=\"text-align:left; valign:top;\">\n"
	

# Begin of form 
	
	data += "<form name=\"date_form\" method=\"POST\" action=\"plot/select\" class=\"command selectfile\">\n"

# Begin of table 
        
	data += "<table border=\"1\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"

# Define cell widths

        data += "<colgroup><col width=\"100\"><col width=\"250\"></colgroup>\n"
	
#     ----------   Date     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\"><b>Date</b> <br> from/to</td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
# Date Year from
        if l_nightly:
	  data += "    <select  disabled=\"disabled\">\n"
	else:
	  data += "    <select>\n"
	  
        data += "      <option style=\"font-size: 8pt\" selected>" + svn_Date_from.split('-')[0] + "</option>\n"
        data += "    </select>\n"

# Date Month from

        if l_nightly:
	  data += "    <select  disabled=\"disabled\">\n"
	else:
	  data += "    <select>\n"

        data += "      <option style=\"font-size: 8pt\" selected>" + svn_Date_from.split('-')[1] + "</option>"
        data += "    </select>\n"
	
# Date Day from
	
	if l_nightly:
	  data += "    <select  disabled=\"disabled\">\n"
	else:
	  data += "    <select>\n"

        data += "      <option  selected>" + svn_Date_from.split('-')[2].split(' ')[0] + "</option>\n"
        data += "    </select>\n"
	
        data += "    <br>\n"

# Date Year to
        if l_nightly:
	  data += "    <select  disabled=\"disabled\">\n"
	else:
	  data += "    <select>\n"

        data += "      <option selected>" + svn_Date_to.split('-')[0] + "</option>"
        data += "    </select>\n"
	
# Date Month to
        if l_nightly:
	  data += "    <select  disabled=\"disabled\">\n"
	else:
	  data += "    <select>\n"

        data += "      <option selected>" + svn_Date_to.split('-')[1] + "</option>"
        data += "    </select>\n"
	
# Date Day from
        if l_nightly:
	  data += "    <select  name=\"Date\" disabled=\"disabled\">\n"
	else:
	  data += "    <select name=\"Dadate_formte\">\n"

        data += "      <option selected>" + svn_Date_to.split('-')[2].split(' ')[0] + "</option>"
        data += "    </select>\n"
	
	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Revision     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\"><b>Rev. Nr.</b><br>from/to</td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
# Revision from 

        if l_nightly:
	  data += "    <select  name=\"rev_from\" disabled=\"disabled\">\n"
	else:
	  data += "    <select  name=\"rev_from\">\n"
	
	data += "    <option selected>" + svn_Rev_from + "</option>\n"
	data += "  </select>\n<br>\n"
        
# Revision to 

        if l_nightly:
	  data += "    <select  name=\"rev_to\" disabled=\"disabled\">\n"
	else:
	  data += "    <select  name=\"rev_to\">\n"
	
	data += "    <option selected>" + svn_Rev_from + "</option>\n"
	data += "  </select>\n"

	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Computer     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Computer</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
        if l_nightly:
	  data += "    <select  name=\"comp\" disabled=\"disabled\">\n"
	else:
	  data += "    <select  name=\"comp\">\n"
	
	data += "    <option selected> all </option>\n"
	data += "  </select>\n"

	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Build Nr.     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Build Nr.</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"

        if l_nightly:
	  data += "    <select  name=\"build_nr\" disabled=\"disabled\">\n"
	else:
	  data += "    <select  name=\"build_nr\">\n"
	
	data += "    <option selected> all </option>\n"
	data += "  </select>\n"

	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Exp. name     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Exp.</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
        if l_nightly:
	  data += "    <select  name=\"exp\" disabled=\"disabled\">\n"
	else:
	  data += "    <select  name=\"exp\">\n"
	
	data += "      <option selected>" + exp_plot_Info + "</option>"
	
	data += "  </select>\n"
	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Exp. file     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>File</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
        data += "    <select  name=\"file\">\n"
	      
        if exp_file_Info == "NotSet" and len(Archive_Button_Dict['file']) > 0:
	  exp_file_Info = Archive_Button_Dict['file'][0].strip(".png")

        for f in Archive_Button_Dict['file']:
	  if f.strip(".png") == exp_file_Info:
	    data += "      <option selected>" + f.strip(".png") + "</option>\n"
	  else:
	    data += "      <option>" + f.strip(".png") + "</option>\n"

	data += "  </select>\n"
	data += "  </td>\n"
	data += "</tr>\n"
        
#     ----------    Search OK Button     ----------   

        data += "<tr>\n  <td style=\"text-align:left;\">\n    <b>Search</b>\n  </td>\n"
        data += "  <td style=\"text-align:left;\">\n"
	data += "    <input type=\"hidden\" name=\"exp\" value=\""+ exp_plot_Info + "\">\n"
	data += "    <input type=\"submit\" name=\"Search_OK\" value=\"OK\" >\n"
        data += "  </td>\n</tr>\n"

# End of ´ICON Buildbot Archive´ and Form area

        data += "</table>\n"
        data += "</form>\n"

# End of ´ICON Buildbot Archive´ div
        
	data += "</div>\n"

#==================================================
#
# Begin of ´Available Plot List´ div
#
#==================================================

        data += "<div id=\"arch_file_name\">\n"
        data += "<h1>Available Plot List</h1>\n"
        data += "<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
	data += "<tr><td style=\"text-align:left; valign:top;\">\n"
	
	e = 0
	r = 0
	d = 0
	c = 0
	b = 0
	f = 0
	

        for Ex in Exp_Dict:
	  e += 1
          data += "\n<!--- " + str(e) + ". Experiment " + Ex + " -->\n"
          data += "\n<a href=\"javascript:anzeigen('date_" + str(e) + "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	  data += "  <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" + Ex + "</span>\n"
	  data += "</a><br>\n"
	  
	  data += "<div id=\"date_" + str(e) + "\" class=\"Rand\" style=\"display: none;\">\n"
	  for Da in Exp_Dict[Ex]:
	    d += 1 
            data += "  <!--- " + str(d) + ". Date " + Da + " -->\n"
          
	    data += "  <a href=\"javascript:anzeigen('rev_" + str(d) + "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	    data += "    <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" + Da + "</span>\n"
	    data += "  </a><br>\n"
#------------------------------
	    data += "  <div id=\"rev_" + str(d) + "\" class=\"Rand\" style=\"display: none;\">\n"
	    for Re in Exp_Dict[Ex][Da]:
	      r += 1 
              data += "    <!--- " + str(r) + ". Revision " + Re + " -->\n"
          
	      data += "    <a href=\"javascript:anzeigen('comp_" + str(r) + "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	      data += "      <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" + Re + "</span>\n"
	      data += "    </a><br>\n"
#------------------------------
	      data += "    <div id=\"comp_" + str(r) + "\" class=\"Rand\" style=\"display: none;\">\n"
	      for Co in Exp_Dict[Ex][Da][Re]:
	        c += 1 
                data += "      <!--- " + str(c) + ". Computer " + Co + " -->\n" 
	        data += "      <a href=\"javascript:anzeigen('build_" + str(c) + "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	        data += "        <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" + Co + "</span>\n"
	        data += "      </a><br>\n"
#------------------------------
	        data += "      <div id=\"build_" + str(c) + "\" class=\"Rand\" style=\"display: none;\">\n"
	        for Bu in Exp_Dict[Ex][Da][Re][Co]:
	          b += 1 
                  data += "        <!--- " + str(b) + ". Build " + Bu + " -->\n" 
	          data += "        <a href=\"javascript:anzeigen('file_" + str(b) + "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	          data += "          <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" + Bu + "</span>\n"
	          data += "        </a><br>\n"
#------------------------------
	          data += "      <div id=\"file_" + str(b) + "\" class=\"Rand\" style=\"display: none;\">\n"
	          for Fi in Exp_Dict[Ex][Da][Re][Co][Bu]:
	            f += 1 
                    data += "          <!--- " + str(f) + ". File " + Co + " -->\n" 
#	            data += "          <a href=\"javascript:anzeigen('file_" + str(f) + "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
#	            data += "            <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" + Fi + "</span>\n"
#	            data += "          </a><br>\n"
                    
		    if Fi.strip(".png") == exp_file_Info:
                      data += "<input type=\"checkbox\" name=\"zutat\" value=\"" + Fi.strip(".png") + "\" id=\"check1\" checked=\"checked\" disabled=\"disabled\">"
                    else:
                      data += "<input type=\"checkbox\" name=\"zutat\" value=\"" + Fi.strip(".png") + "\" id=\"check1\" disabled=\"disabled\">"
		      
		    data += "<label for=\"check1\">"+ Fi.strip(".png") + "</label><br />\n"
#                    data += "<a href=\"#\" class=\"link\"><span style=\"color: #FF0000;\">" + Fi + "</span></a><br>\n"
                  data += "        </div>\n"
#------------------------------
                data += "      </div>\n"
#------------------------------
	      data += "    </div>\n"
#------------------------------
	  
            data += "  </div>\n"
#------------------------------
	  
          data += "</div>\n"
          

        data += "</table>\n"
	  
# End of ´Available Plot List´  div

        data += "</div>\n"

# End of left side  div
        
	data += "</div>\n"
	
#===================================================================================


#===================================================================================

	data += "<div style=\"float:left; padding:3px; margin:5px  width:200px;\">\n"

#===================================================================================

# div vor Reference Plot

#===================================================================================
	
	data += "  <div id=\"arch_ref\">\n"
        
	p = "public_html/reference_plot/" + exp_plot_Info
        if os.path.isdir(p):
	  save_date = os.listdir(p)
	  save_date = sorted(save_date)[-1]
          p += "/" + save_date
	  ref_date = os.listdir(p)[-1]
	  
          p2 = p + "/" + ref_date + "/buildbot/"
	  ref_rev = os.listdir(p2)[-1]

          p3 = p2 + "/" + ref_rev
	  ref_comp = os.listdir(p3)[-1]

          p4 = p3 + "/" + ref_comp
	  ref_build = os.listdir(p4)[-1]
	  
          data += "  <h1>Reference Plot (" + ref_date + ")</h1>\n"
	  data += "  <form name=\"replace_form\" method=\"POST\" action=\"plot/replace\" class=\"command replace\">\n"
          data += "    <input type=\"hidden\" name=\"exp\" value=\""+ exp_plot_Info + "\">\n"
          data += "    <input type=\"submit\" name=\"Replace\" value=\"Replace\" >\n"
	  data += "  </form>\n"
          data += "  <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
	  data += "    <tr>\n    <td>\n"
          data += "      <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
          data += "        <tr><td><b>Rev. Nr:</b></td><td>" + ref_rev + "</td>"
          data += "            <td><b>Computer:</b></td><td>" + ref_comp +"</td>"
          data += "            <td><b>Build. Nr:</b></td><td>" + ref_build + "</td>\n"
          data += "        </tr>\n      </table>\n"
          data += "    </td></tr>\n    <tr><td>\n"
          
	  data += "      <img src=\"reference_plot/"
	  data +=   exp_plot_Info + "/" 
	  data += save_date + "/" 
	  data += ref_date + "/buildbot/" 
	  data += ref_rev + "/"
	  data += ref_comp +"/" 
	  data += ref_build +"/" 
	  data += exp_plot_Info + "/plots/"
	  data += exp_file_Info + ".png\" />\n"
	  
	  data += "    </td></tr>\n"
	  data += "  </table>"
        else:
          data += "<h1>Reference Plot (xxxx-xx-xx)</h1>\n"
	  data += "<form name=\"replace_form\" method=\"POST\" action=\"plot/replace\" class=\"command replace\">"
          data += "  <input type=\"hidden\" name=\"exp\" value=\""+ exp_plot_Info + "\">\n"
          data += "  <input type=\"submit\" name=\"Replace_2\" value=\"Replace\" >\n"
	  data += "</form>"
          data += "<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
	  data += "  <tr>\n    <td>"
          data += "      <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
          data += "        <tr><td><b>Rev. Nr:</b></td><td> ????</td>\n"
          data += "            <td><b>Computer:</b></td><td> ???? </td>\n"
          data += "            <td><b>Build. Nr:</b></td><td> ???? </td>\n"
          data += "</tr>\n      </table>\n"
          data += "    </td></tr>\n    <tr><td>\n"
#	  data += "      <img src=\"archive/ref_" + exp_plot_Info + "/" + ref_date + "/buildbot/" + ref_rev + "/" + ref_comp + "/" + ref_build + "/" + exp_plot_Info + "/plots/" + exp_file_Info + ".png\" />\n"
	  data += "    </td>"
	  data += "  </tr>"
	  data += "</table>"
	  
        print "====== ref Plot ======"
	
#===================================================================================

# div vor Exp Plots

#===================================================================================
	data += "  </div>\n  <div id=\"arch_plots\">\n  <h1>Plot Area</h1>\n"

        p = "public_html/archive/" + svn_Date_from + "/buildbot/" + svn_Rev_from + "/"
        data += "  <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
        
        for comp in Archive_Button_Dict['comp']:
	  for Da in Comp_Dict[comp]:
            for Re in Comp_Dict[comp][Da]:
	      for Bu in Comp_Dict[comp][Da][Re]:
	        for Ex in Comp_Dict[comp][Da][Re][Bu]:
		  for Fi in Comp_Dict[comp][Da][Re][Bu][Ex]:
		    if Fi.strip(".png") == exp_file_Info:
	              data += "  <tr><td>\n"
	              data += "    <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"

	              data += "      <tr>\n"
		      data += "        <td><b>Computer:</b></td><td>" + comp + "</td>\n"
		      data += "        <td><b>Build. Nr:</b></td><td><a href=\"builders/" + comp + "/builds/"+ Bu + "\">" 
		      data += Bu 
		      data += "</a></td></tr>\n"
		
	              data += "    </table>\n"
                      data += "  </td></tr>\n"
                      data += "  <tr><td>\n"        
	              data += "    <img src=\"archive/" + Da + "/buildbot/"+ Re + "/" + comp +"/"+ Bu + "/" + Ex + "/plots/" + Fi + "\" />\n"
	              data += "  </td></tr>\n"
		   
        data += "</table>\n"
        data += "</div>\n"
        data += "</div>\n"
        data += "</div>\n"

#        data += '<hr /><div class="footer">\n'
#        data += '</div>\n'
	data += self.footer(status,request)

        return data

    def select(self, req):
      e = req.args.get("exp",[None])[0]
      f = req.args.get("file",[None])[0]
      return Redirect("../plot?exp="+e+"&file="+f+"&modus=nightly")

    def replace(self, req):
      e = req.args.get("exp",[None])[0]
#      f = req.args.get("file",[None])[0]
      return Redirect("../reference?exp="+e)
    
    def default(self, req):
      return Redirect("../plot")
    
    def getChild(self, path, req):
      if path == "select":
        return self.select(req)
      elif path == "replace":
        return self.replace(req)
      else:
        return self.default(req)

    def get_info(self, request):
        global exp_plot_Info
	global exp_file_Info
        global l_nightly
	
        exp_plot_Info="NotSet"
	exp_file_Info="NotSet"
        l_nightly = True
        
        if "exp" in request.args:
          try:
            exp_plot_Info = request.args["exp"][0]
          except ValueError:
            pass

        if "modus" in request.args:
          try:
	    if request.args["modus"][0] == "archive":
              l_nightly = False
          except ValueError:
            pass
        
	if "file" in request.args:
	  try:
	    exp_file_Info = request.args["file"][0]
	  except ValueError:
	    pass
							     
        return None

    def footer(self, status, req):
        # TODO: this stuff should be generated by a template of some sort
        projectURL = status.getProjectURL()
        projectName = status.getProjectName()
        data = '<hr />\n<div class="footer">\n'

        welcomeurl = self.path_to_root(req) + "index.html"
        data += '  [<a href="%s">welcome</a>]\n' % welcomeurl
        data += "  <br />\n"

        data += '  <a href="http://buildbot.sourceforge.net/">Buildbot</a>'
        data += "-%s " % version
        if projectName:
            data += "working for the "
            if projectURL:
                data += "<a href=\"%s\">%s</a> project.\n" % (projectURL,projectName)
            else:
                data += "%s project.\n" % projectName
        data += "  <br />\n"
        data += ("  Page built: " +
                 time.strftime("%a %d %b %Y %H:%M:%S",
                               time.localtime(util.now()))
                 + "\n")
        data += '</div>\n'

        return data
#================================== NEW ================================
