#================================== NEW ================================

from buildbot.status.web.base import HtmlResource
from twisted.web.util import Redirect
import os, time
from buildbot import version, util
from datetime import datetime, timedelta

class EXP_plot(HtmlResource):
    title = "Exp Plots"
    RevisionNr = "0000"

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
        global modus
        global RevisionNr
	global selected_builder
	global selected_build

	global selected_Date_from
	global selected_Date_to
	
	self.get_info(request)
#ToDo when archive has no plots
	if exp_plot_Info == "test_couple_atmo-ocean":
	  data = "<h1>For experiment <b>" + exp_plot_Info + "</b> are no plots available"
	  return data
	
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
#	      if e == exp_plot_Info or exp_plot_Info == "all":
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
            if (r >= Rev_from) and (r <= Rev_to):
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
	         if r in Exp_Dict[e][d]:
	           if c in Exp_Dict[e][d][r]:
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
             if d in Comp_Dict[c]:
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

	RevisionNr = "0000"
        file = open("icon_run/svn_info")
        while 1:
          line = file.readline()
          if not line:
            break
          
          if line.find("Revision:") >= 0:
            tmp,tmp_Rev = line.split("Revision:",1)
            RevisionNr = tmp_Rev.strip(" ").replace("\n","")
            
        file.close

        if l_nightly:
          Rev_from = RevisionNr
	  Rev_to = Rev_from
 
	  Date_to = time.strftime("%Y-%m-%d",time.localtime(util.now()))
#ws	  Date_to = "2012-06-21"
	  now_time = time.strftime("%H%M",time.localtime(util.now()))
	      
	  yesterday = datetime.now() - timedelta(days=2)
	      
	  if now_time > "2200":
#ws	    Date_from = "2012-06-20"
            Date_from = Date_to
	  else:
            Date_from = yesterday.strftime("%Y-%m-%d")
#ws	    Date_from = "2012-06-20"
	
	else:
          Rev_from  = str(int(RevisionNr)-50)
          Rev_to    = RevisionNr
          now_min_10d = datetime.now() - timedelta(days=10)
          Date_from = now_min_10d.strftime("%Y-%m-%d")
          Date_to = time.strftime("%Y-%m-%d",time.localtime(util.now()))
          
          if selected_Date_from != "NotSet":
	    Date_from = selected_Date_from  
	    
          if selected_Date_to != "NotSet":
	    Date_to = selected_Date_to

	  if selected_Rev_from != "NotSet":
	    Rev_from = selected_Rev_from  
	    
          if selected_Rev_to != "NotSet":
	    Rev_to = selected_Rev_to

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
	  if DATE.find("20") == 0:
            if (DATE >= Date_from) and (DATE <= Date_to):
	      if add_rev(p_date,DATE,date_Dict):
	        Archive_Button_Dict['date'].append(DATE)
	    
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

        Archive_Button_Dict['date']  = list(sorted(set(Archive_Button_Dict['date'])))
 	Archive_Button_Dict['rev']   = list(sorted(set(Archive_Button_Dict['rev'])))
	Archive_Button_Dict['comp']  = list(sorted(set(Archive_Button_Dict['comp'])))
	Archive_Button_Dict['build'] = list(sorted(set(Archive_Button_Dict['build'])))
	Archive_Button_Dict['exp']   = list(sorted(set(Archive_Button_Dict['exp'])))
	Archive_Button_Dict['file']  = list(sorted(set(Archive_Button_Dict['file'])))
	
        if not Archive_Button_Dict['rev']:
	  data = "Revision List for this timeperode ist emty"
	  return data
	  
	print Archive_Button_Dict['rev']
        if not l_nightly:
	  min_Revision = Archive_Button_Dict['rev'][0]
	  max_Revision = Archive_Button_Dict['rev'][-1]
	  
	  if Rev_from < min_Revision:
	    Rev_from = min_Revision
	    
	  if Rev_to > max_Revision:
	    Rev_to   = max_Revision
	    
	

# Build a plot list which contains only the plots of the acual experiment

	t_file = []
	for Fi in Archive_Button_Dict['file']:
	  if Fi.find(exp_plot_Info+"_") >= 0:
	    t_file.append(Fi) 
	    print Fi + " (" + exp_plot_Info + "_) included in the list"
	  else:
	    print Fi + " (" + exp_plot_Info + "_) not included in the list"
	    
        Archive_Button_Dict['file'] = t_file
	
#ToDo	if len(Archive_Button_Dict['comp']) > 1:
#ToDo	  Archive_Button_Dict['comp'].insert(0, 'all')

#========================================================================================================
#
# End searching for plots
#
#========================================================================================================

#==============================================
# Begin of whole div
        

        data = ""
        
	
#ToDo        if l_nightly:
#ToDo	  data +=  "<table>\n"
#ToDo	  data +=  "  <tr>\n"
#ToDo	  data +=  "    <td id=\"exp_header\"> Source:</td><td id=\"exp_header\"> trunk/icon-dev &#64; r" + RevisionNr  + "</td>\n"
#ToDo	  data +=  "  </tr><tr>\n"
#ToDo	  data +=  "    <td id=\"exp_header\">Experiment:</td><td id=\"exp_header\">" + exp_plot_Info + "</td>\n"
#ToDo	  data +=  "  </tr><tr>\n"
#ToDo	  data +=  "    <td id=\"exp_header\">Experiment descriptor file: </td>"
#ToDo	  data +=  "    <td id=\"exp_header\"><a href=\""
#ToDo	  data +=  " https://code.zmaw.de/projects/icon/repository/entry/trunk/icon-dev/run/exp." + exp_plot_Info 
#ToDo	  data += "?rev="+RevisionNr+"\">exp." + exp_plot_Info + "</a></td>\n"
#ToDo	  data +=  "  </tr><tr>\n"
#ToDo	  data +=  "    <td id=\"exp_header\">Post-processing descriptor file: </td>"
#ToDo	  data +=  "    <td id=\"exp_header\"><a href=\""
#ToDo	  data +=  " https://code.zmaw.de/projects/icon/repository/entry/trunk/icon-dev/run/post." + exp_plot_Info 
#ToDo	  data +=  "?rev="+RevisionNr + "\">post." + exp_plot_Info + "</a></td>\n"
#ToDo	  data +=  "  </tr>\n"
#ToDo	  data +=  "</table>\n"
#ToDo     data +=  "<br />\n"
          
        data +=  "<div id=\"page_menu\" >\n"

# Begin of left side  div
	
	data += "<div style=\"float:left; padding:3px; margin:5px width:350px;\">\n"

# Begin of ´ICON Buildbot Archive´ div
        
        data += "<div id=\"arch_menu\" >\n"
	data += "<h1>ICON Buildbot Archive</h1>\n"
	

# Begin of form 
	
	data += "<form name=\"date_form\" method=\"POST\" action=\"plot/select\" class=\"command selectfile\">\n"

# Begin of table 
        
	data += "<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"540\">\n"

# Define cell widths

        data += "<colgroup><col width=\"80\"><col width=\"470\"></colgroup>\n"
	
#     ----------   Date     ----------   

        data += "<tr>\n"
#	data += "  <td style=\"text-align:left;\"><b>Date:</b> <br> from/to</td>\n"
	data += "  <td style=\"text-align:left;\"><b>Date:</b></td>\n"
	data += "  <td style=\"text-align:left;\">"
	
# Date Year from
        if l_nightly:
	  data += Date_from 
	else:
	  data += "    <input name=\"date_from\" type=\"text\" value=\"" + Date_from.split(' ')[0] + "\" size=\"10\" maxlength=\"10\" >\n"
	  

# Date Year to
        if l_nightly:
	  if Date_to != Date_from:
	    data += "&nbsp; to &nbsp" + Date_to
	else:
	  data += "&nbsp; to &nbsp"
	  data += "    <input name=\"date_to\" type=\"text\" value=\"" + Date_to.split(' ')[0] + "\" size=\"10\" maxlength=\"10\">\n"
	  data += "    <input type=\"submit\" name=\"Date_button\" value=\">\">"
	
	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Revision     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\"><b>Revision:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
# Revision from 

        if l_nightly:
	  data += Rev_from
	else:
	  data += "    <input name=\"rev_from\" type=\"text\" value=\"" + Rev_from.split(' ')[0] + "\" size=\"10\" maxlength=\"10\">\n"
        
# Revision to 

        if l_nightly:
	  if Rev_to != Rev_from:
	    data += "&nbsp; to &nbsp" + Rev_to
	else:
	  data += "&nbsp; to &nbsp"
	  data += "    <input name=\"rev_to\" type=\"text\" value=\"" + Rev_to.split(' ')[0] + "\" size=\"10\" maxlength=\"10\">\n"
	  data += "    <input type=\"submit\" name=\"Rev_button\" value=\">\">"

	data += "  </td>\n"
	data += "</tr>\n"
	
# Branch 

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>SVN-Branch:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"

        if l_nightly:
	  data += "trunk/icon-dev"
	else:
	  data += "    <select  name=\"branch\" onChange=\"document.date_form.submit()\">\n"
	  if selected_build == "all":
	    data += "    <option selected> all </option>\n"
	  else:
	    data += "    <option> all </option>\n"
	    
#	  for Bran in Archive_Button_Dict['branch']:
#	    if selected_branch == Bu:
#	      data += "    <option selected>"+ Bran + "</option>\n"
#	    else:
#	      data += "    <option>"+ Bran + "</option>\n"
	      
	  data += "  </select>\n"

	data += "  </td>\n"
	data += "</tr>\n"

#     ----------    Computer     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Builder:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
        if l_nightly:
	  data += "all"
	else:
	  data += "    <select  name=\"builder\" onChange=\"document.date_form.submit()\">\n"
	  if selected_builder == "all":
	    data += "    <option selected> all </option>\n"
	  else:
	    data += "    <option> all </option>\n"
	
	  for Co in Archive_Button_Dict['comp']:
	    if selected_builder == Co:
	      data += "    <option selected>"+ Co + "</option>\n"
	    else:
	      data += "    <option>"+ Co + "</option>\n"
	  data += "  </select>\n"

	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Build Nr.     ----------   
	t_build = []
        if selected_builder != "all":
	  for Da in Comp_Dict[selected_builder]:
            for Re in Comp_Dict[selected_builder][Da]:
	      for Bu in Comp_Dict[selected_builder][Da][Re]:
		t_build.append(Bu)
	
        Archive_Button_Dict['build'] = list(sorted(t_build))

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Build:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"

        if l_nightly:
	  data += "all"
	else:
	  data += "    <select  name=\"build\" onChange=\"document.date_form.submit()\">\n"
	  if selected_build == "all":
	    data += "    <option selected> all </option>\n"
	  else:
	    data += "    <option> all </option>\n"
	    
	  for Bu in Archive_Button_Dict['build']:
	    if selected_build == Bu:
	      data += "    <option selected>"+ Bu + "</option>\n"
	    else:
	      data += "    <option>"+ Bu + "</option>\n"
	      
	  data += "  </select>\n"

	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Exp. name     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Experiment:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
	data += "    <select  name=\"exp\" onChange=\"document.date_form.submit()\">\n"
#	data += "    <select  name=\"exp\">\n"

#ToDo   Include the correct showing of '_' in experment name and plot-name

	for Ex in Archive_Button_Dict['exp']:
	  if Ex == exp_plot_Info:
	    data += "      <option selected>" + Ex + "</option>\n"
	  else:
	    data += "      <option>" + Ex + "</option>\n"
	
	
	data += "  </select>\n"
	data += "  </td>\n"
	data += "</tr>\n"
	
#     ----------    Exp. file     ----------   

        data += "<tr>\n"
	data += "  <td style=\"text-align:left;\" ><b>Plot:</b></td>\n"
	data += "  <td style=\"text-align:left;\">\n"
	
        data += "    <select  name=\"file\" onChange=\"document.date_form.submit()\">\n"
	      
        if exp_file_Info == "NotSet" and len(Archive_Button_Dict['file']) > 0:
	  exp_file_Info = Archive_Button_Dict['file'][0].strip(".png")

        for f in Archive_Button_Dict['file']:
	  if f.strip(".png") == exp_file_Info:
	    data += "      <option selected>" + f.strip(".png") + "</option>\n"
	  else:
	    data += "      <option>" + f.strip(".png") + "</option>\n"

	data += "    </select>\n"
	data += "  </td>\n"
	data += "</tr>\n"
        
#     ----------    Search OK Button     ----------   

#        data += "<tr>\n  <td style=\"text-align:left;\">\n    <b>Search</b>\n  </td>\n"
#        data += "  <td style=\"text-align:left;\">\n"
	data += "    <input type=\"hidden\" name=\"modus\" value=\""+ modus + "\">\n"
	data += "    <input type=\"hidden\" name=\"exp\" value=\""+ exp_plot_Info + "\">\n"
	data += "    <input type=\"hidden\" name=\"builder\" value=\"" + selected_builder + "\">\n"
	data += "    <input type=\"hidden\" name=\"build\" value=\"" + selected_build + "\">\n"
#      fDate   = req.args.get("date_from",[None])[0]
#      tDate   = req.args.get("date_to",[None])[0]
#      fRev    = req.args.get("rev_from",[None])[0]
#      tRev    = req.args.get("rev_to",[None])[0]
#	data += "    <input type=\"submit\" name=\"Search_OK\" value=\"OK\" >\n"
#        data += "  </td>\n</tr>\n"

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
	
        for Ex in Archive_Button_Dict['exp']:
#        for Ex in Exp_Dict:
	  e += 1
          data += "\n<!--- " + str(e) + ". Experiment " + Ex + " -->\n"
          data += "\n<a href=\"javascript:anzeigen('date_" + str(e) 
	  data += "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	  data += "  <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" + Ex + "</span>\n"
	  data += "</a><br>\n"
	  
	  data += "<div id=\"date_" + str(e) + "\" class=\"Rand\" style=\"display: none;\">\n"
	  for Da in Exp_Dict[Ex]:
	    d += 1 
            data += "  <!--- " + str(d) + ". Date " + Da + " -->\n"
          
	    data += "  <a href=\"javascript:anzeigen('rev_" + str(d) 
	    data += "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	    data += "    <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" 
	    data += Da + "</span>\n"
	    data += "  </a><br>\n"
#------------------------------
	    data += "  <div id=\"rev_" + str(d) + "\" class=\"Rand\" style=\"display: none;\">\n"
	    for Re in Exp_Dict[Ex][Da]:
	      r += 1 
              data += "    <!--- " + str(r) + ". Revision " + Re + " -->\n"
          
	      data += "    <a href=\"javascript:anzeigen('comp_" + str(r) 
	      data += "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	      data += "      <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" 
	      data += Re + "</span>\n"
	      data += "    </a><br>\n"
#------------------------------
	      data += "    <div id=\"comp_" + str(r) + "\" class=\"Rand\" style=\"display: none;\">\n"
	      for Co in Exp_Dict[Ex][Da][Re]:
	        c += 1 
                data += "      <!--- " + str(c) + ". Computer " + Co + " -->\n" 
	        data += "      <a href=\"javascript:anzeigen('build_" + str(c) 
		data += "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	        data += "        <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" 
		data += Co + "</span>\n"
	        data += "      </a><br>\n"
#------------------------------
	        data += "      <div id=\"build_" + str(c) + "\" class=\"Rand\" style=\"display: none;\">\n"
	        for Bu in Exp_Dict[Ex][Da][Re][Co]:
	          b += 1 
                  data += "        <!--- " + str(b) + ". Build " + Bu + " -->\n" 
	          data += "        <a href=\"javascript:anzeigen('file_" + str(b) 
		  data += "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	          data += "          <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" 
		  data +=  Bu + "</span>\n"
	          data += "        </a><br>\n"
#------------------------------
	          data += "      <div id=\"file_" + str(b) + "\" class=\"Rand\" style=\"display: none;\">\n"
	          for Fi in Exp_Dict[Ex][Da][Re][Co][Bu]:
	            f += 1 
                    data += "          <!--- " + str(f) + ". File " + Co + " -->\n" 
                    
                    data += "<input type=\"checkbox\" name=\"plot_check\" value=\"" + Fi.strip(".png") + "\" id=\"check1\""
		    if Fi.strip(".png") == exp_file_Info:
		      data += "checked=\"checked\""
		      
		    if l_nightly:
		      data += "disabled=\"disabled\""
		    data += ">"
		    data += "<label for=\"check1\">"+ Fi.strip(".png") + "</label><br />\n"
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
	
        if l_nightly:
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
	  
            data += "  <h1>Reference Plot </h1>\n"
            
	    data += "  <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"600\">\n"
	    data += "    <tr>\n    <td>\n"
            data += "      <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"600\">\n"
            data += "        <colgroup>"
	    data += "          <col width=\"60\">"
	    data += "          <col width=\"40\">"
	    data += "          <col width=\"60\">"
	    data += "          <col width=\"120\">"
	    data += "          <col width=\"50\">"
	    data += "          <col width=\"50\">"
	    data += "          <col width=\"60\">"
	    data += "          <col width=\"70\">"
	    data += "          <col width=\"90\">"
	    data += "        </colgroup>\n"
            data += "        <tr><td><b>Revision:</b></td><td>" + ref_rev + "</td>"
            data += "            <td><b>Builder:</b></td><td>" + ref_comp +"</td>"
#            data += "            <td><b>Build:</b></td><td>" + ref_build + "</td>"
            data += "        <td><b>Build:</b></td><td><a href=\"builders/" + ref_comp + "/builds/"+ ref_build + "\">" 
	    data += ref_build 
	    data += "</a></td>\n"
            data += "            <td><b>Date:</b></td><td>" + ref_date + "</td>\n"
            
	    data += "            <td style=\"text-align:left;\">\n"
	    data += "  <form name=\"replace_form\" method=\"POST\" action=\"plot/replace\" class=\"command replace\">\n"
            data += "    <input type=\"hidden\" name=\"exp\" value=\""+ exp_plot_Info + "\">\n"
            data += "    <input type=\"hidden\" name=\"rev\" value=\""+ ref_rev + "\">\n"
            data += "    <input type=\"submit\" name=\"Replace\" value=\"Replace\" >\n"
	    data += "  </form>\n"
            
	    data += "            </td>\n"

	    data += "        </tr>\n      </table>\n"
            data += "    </td></tr>\n    <tr><td style=\"text-align:left; valign:top;\">\n"
          
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
            data += "<h1>Reference Plot</h1>\n"
	    data += "<form name=\"replace_form\" method=\"POST\" action=\"plot/replace\" class=\"command replace\">"
            data += "  <input type=\"hidden\" name=\"exp\" value=\""+ exp_plot_Info + "\">\n"
            data += "  <input type=\"submit\" name=\"Replace_2\" value=\"Replace\" >\n"
	    data += "</form>"
            data += "<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"600\">\n"
	    data += "  <tr>\n    <td>"
            data += "      <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"600\">\n"
            data += "        <colgroup>"
	    data += "          <col width=\"60\">"
	    data += "          <col width=\"40\">"
	    data += "          <col width=\"60\">"
	    data += "          <col width=\"120\">"
	    data += "          <col width=\"50\">"
	    data += "          <col width=\"50\">"
	    data += "          <col width=\"60\">"
	    data += "          <col width=\"70\">"
	    data += "          <col width=\"90\">"
	    data += "        </colgroup>\n"
            data += "        <tr><td><b>Revision:</b></td><td> ????</td>\n"
            data += "            <td><b>Builder:</b></td><td> ???? </td>\n"
            data += "            <td><b>Build:</b></td><td> ???? </td>\n"
            data += "            <td><b>Date:</b></td><td> xxxx-xx-xx </td>\n"
            data += "</tr>\n      </table>\n"
            data += "    </td></tr>\n    <tr><td style=\"text-align:left; valign:top;\">\n"
	    data += "    </td>"
	    data += "  </tr>"
	    data += "</table>"
	    
	  data += "</div>\n "
	  
          print "====== ref Plot ======"
	
#===================================================================================

# div vor Exp Plots

#===================================================================================
        if l_nightly:
	  data += "<div id=\"ap_ref\">\n  <h1>Plot Area</h1>\n"
	else:
	  data += "<div id=\"arch_plots\">\n  <h1>Plot Area</h1>\n"
	

        p = "public_html/archive/" + Date_from + "/buildbot/" + Rev_from + "/"
        data += "  <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"600\">\n"
        
        for comp in Archive_Button_Dict['comp']:
	  if comp == selected_builder or selected_builder == "all":
	    for Da in sorted(Comp_Dict[comp]):
              for Re in sorted(Comp_Dict[comp][Da]):
	        for Bu in sorted(Comp_Dict[comp][Da][Re]):
	          if Bu == selected_build or selected_build == "all":
	            for Ex in Comp_Dict[comp][Da][Re][Bu]:
		      for Fi in Comp_Dict[comp][Da][Re][Bu][Ex]:
		        if Fi.strip(".png") == exp_file_Info:
	                  data += "  <tr><td>\n"
	                  data += "    <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"600\">\n"
            
	                  data += "        <colgroup>"
	                  data += "          <col width=\"60\">"
	                  data += "          <col width=\"40\">"
	                  data += "          <col width=\"60\">"
	                  data += "          <col width=\"120\">"
	                  data += "          <col width=\"50\">"
	                  data += "          <col width=\"50\">"
	                  data += "          <col width=\"60\">"
	                  data += "          <col width=\"70\">"
	                  data += "          <col width=\"90\">"
	                  data += "        </colgroup>\n"
               
	                  data += "      <tr>\n"
                          data += "        <td><b>Revision:</b></td><td>" + Re + "</td>"
		          data += "        <td><b>Builder:</b></td><td>" + comp + "</td>\n"
		          data += "        <td><b>Build:</b></td><td><a href=\"builders/" + comp + "/builds/"+ Bu + "\">" 
		          data += Bu 
		          data += "</a></td>\n"
                          data += "        <td><b>Date:</b></td><td>" + Da + "</td>\n"
                          data += "        <td><b>&nbsp</b></td><td>&nbsp</td>\n"
		          data += "</tr>\n"
		
	                  data += "    </table>\n"
                          data += "  </td></tr>\n"
                          data += "  <tr><td style=\"text-align:left;\">\n"        
	                  data += "    <img src=\"archive/" + Da + "/buildbot/"+ Re + "/" + comp +"/"
		          data +=  Bu + "/" + Ex + "/plots/" + Fi + "\" />\n"
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
      global selected_builder
      print "======  select ====="
      print req
      print req.args
      print "======  select ====="
      m       = req.args.get("modus",[None])[0]
      e       = req.args.get("exp",[None])[0]
      f       = req.args.get("file",[None])[0]
      builder = req.args.get("builder",[None])[0]
      build   = req.args.get("build",[None])[0]
#      fDate   = req.args.get("date_from",[None])[0]
#      tDate   = req.args.get("date_to",[None])[0]
#      fRev    = req.args.get("rev_from",[None])[0]
#      tRev    = req.args.get("rev_to",[None])[0]
      
      Redirect_Info = "../plot?exp=" + e  + "&modus=" + m
      
      if f.find(e + "_") >= 0:
        Redirect_Info += "&file=" + f

      Redirect_Info += "&build=" + build
      
      if builder != "all":
	Redirect_Info += "&builder=" + builder
	if selected_builder != builder:
          build = "all"
          
      if "date_from" in req.args:
        try:
          fDate = req.args["date_from"][0]
          Redirect_Info += "&date_from=" + fDate      
        except ValueError:
          pass
      
      if "date_to" in req.args:
        try:
          tDate = req.args["date_to"][0]
          Redirect_Info += "&date_to=" + tDate      
        except ValueError:
          pass
      
      if "rev_from" in req.args:
        try:
          fRev = req.args["rev_from"][0]
          Redirect_Info += "&rev_from=" + fRev      
        except ValueError:
          pass
      
      if "rev_to" in req.args:
        try:
          tRev = req.args["rev_to"][0]
          Redirect_Info += "&rev_to=" + tRev      
        except ValueError:
          pass

      return Redirect(Redirect_Info)
#      return Redirect("../plot?exp=" + e + "&file=" + f + "&modus=")

    def replace(self, req):
      e = req.args.get("exp",[None])[0]
      r = req.args.get("rev",[None])[0]
#      f = req.args.get("file",[None])[0]
      return Redirect("../reference?exp="+e+"&start="+r)
    
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
        global modus
	global selected_builder
	global selected_build

	global selected_Date_from
	global selected_Date_to
	
	global selected_Rev_from
	global selected_Rev_to

	exp_plot_Info      = "NotSet"
	exp_file_Info      = "NotSet"
	modus              = "nightly"
	selected_builder   = "all"
	selected_build     = "all"
	selected_Date_from = "NotSet"
	selected_Date_to   = "NotSet"
	selected_Rev_from  = "NotSet"
	selected_Rev_to    = "NotSet"
	
        l_nightly = True
        
        if "exp" in request.args:
          try:
            exp_plot_Info = request.args["exp"][0]
          except ValueError:
            pass

        if "builder" in request.args:
          try:
            selected_builder = request.args["builder"][0]
          except ValueError:
            pass
        
        if "build" in request.args:
          try:
            selected_build = request.args["build"][0]
          except ValueError:
            pass

        if "modus" in request.args:
          try:
	    if request.args["modus"][0] == "archive":
	      modus = "archive"
              l_nightly = False
          except ValueError:
            pass
        
	if "file" in request.args:
	  try:
	    exp_file_Info = request.args["file"][0]
	  except ValueError:
	    pass
							     
	if "date_from" in request.args:
	  try:
	    selected_Date_from = request.args["date_from"][0]
	  except ValueError:
	    pass

	if "date_to" in request.args:
	  try:
	    selected_Date_to = request.args["date_to"][0]
	  except ValueError:
	    pass

	if "rev_from" in request.args:
	  try:
	    selected_Rev_from = request.args["rev_from"][0]
	  except ValueError:
	    pass

	if "rev_to" in request.args:
	  try:
	    selected_Rev_to = request.args["rev_to"][0]
	  except ValueError:
	    pass

	return None

    def footer(self, status, req):
        global RevisionNr
        # TODO: this stuff should be generated by a template of some sort
        projectURL = status.getProjectURL()
        projectName = status.getProjectName()
        data = '<div class="footer">\n'

        welcomeurl = self.path_to_root(req) + "index.html"
        data += '  [<a href="%s">mainpage</a>]\n' % welcomeurl
        data += "  <br />\n"

#        data += "  <br />\n"
        data += ("  Page built: " +
                 time.strftime("%a %d %b %Y %H:%M:%S",
                               time.localtime(util.now()))
                 + "\n")
        data += '</div>\n'

        return data
#================================== NEW ================================
