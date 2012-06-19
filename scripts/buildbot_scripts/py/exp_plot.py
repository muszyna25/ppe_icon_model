#================================== NEW ================================

from buildbot.status.web.base import HtmlResource
import os

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
        global l_nightly
	self.get_info(request)
	
	Archive_Button_Dict = {}
	Archive_Button_Dict['date']  = []
	Archive_Button_Dict['rev']   = []
	Archive_Button_Dict['comp']  = []
	Archive_Button_Dict['build'] = []
	Archive_Button_Dict['exp']   = []
	Archive_Button_Dict['file']  = []
	
#        ExpListInfo = []

	
        file = open("icon_run/svn_info")
        while 1:
          line = file.readline()
          if not line:
            break
          
          if line.find("Revision:") >= 0:
	    tmp,tmp_Rev = line.split("Revision:",1)
            svn_Rev = tmp_Rev.strip(" ").replace("\n","")
                
          if line.find("Last Changed Date:") >= 0:
	    tmp,tmp_Date = line.split("Last Changed Date:",1)
            tmp_Date.replace("\n","")
            tmp_Date.strip(" ")
            svn_Date = tmp_Date.split(' ')[1].replace("\n","")
               
        file.close
         
         
# Build new emty date_Dict type
        date_Dict = {}

# Build new emty rev_Dict type and save it in the date_Dict under the correct Date
        rev_Dict = {}
        date_Dict[svn_Date] = rev_Dict
	Archive_Button_Dict['date'].append(svn_Date)

# Build new emty comp_Dict type and save it in the rev_Dict under the correct reverents number
        comp_Dict = {}
	rev_Dict[svn_Rev] = comp_Dict
	Archive_Button_Dict['rev'].append(svn_Rev)
	
        p = "./public_html/archive/" + svn_Date + "/buildbot/" + svn_Rev + "/"
	for computer in os.listdir(p):
# Build for each computer a new emty build_Dict type and save it in the comp_Dict under the correct computer name
          build_Dict =  {}
          comp_Dict[computer] = build_Dict
  	  Archive_Button_Dict['comp'].append(computer)
          for build in os.listdir(p + "/" + computer + "/"):
# Build for each build number a new emty exp_Dict type and save it in the build_Dict under the correct builder number
            exp_Dict = {}
            build_Dict[build] = exp_Dict
    	    Archive_Button_Dict['build'].append(build)
            for exp in os.listdir(p + "/" + computer + "/" + build + "/"):
	      if exp != "run_info.txt": 
# To Do
# Hier muss noch eingebaut werden das bei 'all' alle Experimente genommen werden.
		if exp == exp_plot_Info:
	          file_List = []
                  exp_Dict[exp] = file_List
      	          Archive_Button_Dict['exp'].append(exp)
# Build a new file list and include the all stateman into the list
                  for fi in os.listdir(p + "/" + computer + "/" + build + "/" + exp + "/plots/"):
# Append all existing files from the the experiment
	            file_List.append(fi)
       	            Archive_Button_Dict['file'].append(fi)

        Archive_Button_Dict['date']  = list(sorted(set(Archive_Button_Dict['date'])))
 	Archive_Button_Dict['rev']   = list(sorted(set(Archive_Button_Dict['rev'])))
	Archive_Button_Dict['comp']  = list(sorted(set(Archive_Button_Dict['comp'])))
	Archive_Button_Dict['build'] = list(sorted(set(Archive_Button_Dict['build'])))
	Archive_Button_Dict['exp']   = list(sorted(set(Archive_Button_Dict['exp'])))
	Archive_Button_Dict['file']  = list(sorted(set(Archive_Button_Dict['file'])))
	
	if len(Archive_Button_Dict['comp']) > 1:
	  Archive_Button_Dict['comp'].insert(0, 'all')
#	if len(Archive_Button_Dict['build']) > 1:
#	  Archive_Button_Dict['build'].insert(0, 'all')
#	if len(Archive_Button_Dict['exp']) > 1:
#	  Archive_Button_Dict['exp'].insert(0, 'all')
#	if len(Archive_Button_Dict['file']) > 1:
#	  Archive_Button_Dict['file'].insert(0, 'all')
	
        data = '''
  <div style="float:left; padding:3px; margin:5px width:350px;">
    <div id="arch_menu" >
      <h1>ICON Buildbot Archive</h1>
      <table border="0" cellspacing="0" cellpadding="0" align="left" width="400">
        <colgroup>
          <col width="100">
          <col width="250">
        </colgroup>

        <tr>
          <td><b>Date</b> <br> from/to</td>
          <td>
	    <form name=date_form>
	    
	'''
        if l_nightly == "nightly":
	  data += "<select  disabled=\"disabled\">\n"
	else:
	  data += "<select>\n"
	  
        data += "<option selected>" + tmp_Date.split('-')[0] + "</option></select>\n"
        data += "</select>\n"
	
        if l_nightly == "nightly":
	  data += "<select  disabled=\"disabled\">\n"
	else:
	  data += "<select>\n"

        data += "<option selected>" + tmp_Date.split('-')[1] + "</option></select>"
        data += "</select>\n"
	
	data += "<select  disabled=\"disabled\">\n"
        data += "<option selected>" + tmp_Date.split('-')[2].split(' ')[0] + "</option></select>\n"
        data += "</select>\n"
	
	data += '''
	    </form>
	    <form name=date_form>
	'''

        if l_nightly == "nightly":
	  data += "<select  disabled=\"disabled\">\n"
	else:
	  data += "<select>\n"

        data += "<option selected>" + tmp_Date.split('-')[0] + "</option></select>"
        data += "</select>\n"
	
        if l_nightly == "nightly":
	  data += "<select  disabled=\"disabled\">\n"
	else:
	  data += "<select>\n"

        data += "<option selected>" + tmp_Date.split('-')[1] + "</option></select>"
        data += "</select>\n"
	
        if l_nightly == "nightly":
	  data += "<select  disabled=\"disabled\">\n"
	else:
	  data += "<select>\n"

        data += "<option selected>" + tmp_Date.split('-')[2].split(' ')[0] + "</option></select>"
        data += "</select>\n"
	
	data += '''
	    </form>
          </td>
        </tr>
        <tr>
          <td><b>Rev. Nr.</b><br>from/to</td>
            <td>
	      <form name=ref_form>
		<select name=mytextarea disabled="disabled">
	'''
	data += "<option selected>" + svn_Rev + "</option></select></form>"
	
	data += "<form name=ref_form >"
	data += "  <select  disabled=\"disabled\">\n"
	data += "    <option selected>" + svn_Rev + "</option>\n"
	data += "  </select>\n"
	data += "</form>"
	
	data += '''
            </td>
          </tr>
          <tr>
            <td width="100px"><b>Computer</b></td>
            <td>
	      <form name=comp_form>
		<select name=mytextarea disabled="disabled">
		  <option selected> all </option>
		</select>
	    </form>
          </td>
        </tr>
        <tr>
          <td><b>Build Nr.</b></td>
          <td>
	    <form name=build_form>
	      <select name=mytextarea disabled="disabled">
		<option  selected> all </option>
	      </select>
	    </form>
          </td>
        </tr>
        <tr>
          <td><b>Exp.</b></td>
          <td>
	    <form name=exp_form >
	      <select name=mytextarea disabled="disabled">
	      '''
	      
	data += "<option selected>" + exp_plot_Info + "</option>"
	
	data += '''
	      </select>
	    </form>
          </td>
        </tr>
        <tr>
          <td><b>File</b></td>
            <td>
	      <form name=file_form>
		<select>
	      '''
        exp_file_Info = Archive_Button_Dict['file'][0]

        for f in Archive_Button_Dict['file']:
	  if f.strip(".png") == exp_file_Info:
	    data += "<option selected>" + f.strip(".png") + "</option>"
	  else:
	    data += "<option>" + f.strip(".png") + "</option>"

	data += '''
	        </select>
	      </form>
            </td>
          </td>
        </tr>
	'''
        
        data += "<tr><td><b>Search</b><td>\n"
        data += "  <input type=\"button\" name=\"Text 1\" value=\"OK\" >\n"
        data += "</td></tr>\n"
        
	data += '''
      </table>
    </div>
    <div id="arch_file_name">
      <h1>Available Plot List</h1>
      <table border="0" cellspacing="0" cellpadding="0" align="left" width="400">
        <tr>
          <td valign=top>

	'''
	e = 0
	r = 0
	d = 0
	c = 0
	b = 0
	f = 0
	
        p = "./public_html/archive/" + svn_Date + "/buildbot/" + svn_Rev + "/"
        e += 1

        r += 1
        data += "\n<!--- " + str(e) + ". Experiment -->\n"

        data += "\n<a href=\"javascript:anzeigen('rev_" + str(r) + "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
	data += "  <img name=\"bild_1\" src=\"closed.gif\" class=\"Bild\"> <span style=\"color: #000000;\">" + exp_plot_Info + "</span>\n"
	data += "</a><br>\n"
	
	data += "<div id=\"rev_" + str(r) + "\" class=\"Rand\" style=\"display: none;\">\n"
        data += "\n<!--- " + str(r) + ". Revision Nr. -->\n"

        d += 1
        data += "  <a href=\"javascript:anzeigen('date_" + str(d) + "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
        data += "    <img name=\"bild_1\" src=\"closed.gif\" lass=\"Bild\"> <span style=\"color: #000000;\">" + svn_Rev + "</span>\n"
        data += "  </a><br>\n"

        data += "  <div id=\"date_" + str(d) + "\" class=\"Rand\" style=\"display: none;\">\n"
        data += "<!--- " + str(d) + ". Datum. -->\n"
        
        c += 1
        data += "    <a href=\"javascript:anzeigen('comp_" + str(c) + "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
        data += "      <img name=\"bild_1\" src=\"closed.gif\" lass=\"Bild\"> <span style=\"color: #000000;\">" + svn_Date + "</span>\n"
        data += "    </a><br>\n"
        data += "    <div id=\"comp_" + str(c) + "\" class=\"Rand\" style=\"display: none;\">\n"
        data += "<!--- " + str(c) + ". Computer -->\n"
        
        for computer in os.listdir(p):
	  b += 1
          data += "    <a href=\"javascript:anzeigen('build_" + str(b) + "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
          data += "      <img name=\"bild_1\" src=\"closed.gif\" lass=\"Bild\"> <span style=\"color: #000000;\">" + computer + "</span>\n"
          data += "    </a><br>\n"
          data += "    <div id=\"build_" + str(b) + "\" class=\"Rand\" style=\"display: none;\">\n"
          data += "<!--- " + str(b) + ". build -->\n"

          for build in os.listdir(p + "/" + computer + "/"):
  	    f += 1
            data += "      <a href=\"javascript:anzeigen('file_" + str(f) + "','bild_1');\" class=\"Ordner\" style=\"text-decoration: None;\">\n"
            data += "        <img name=\"bild_1\" src=\"closed.gif\" lass=\"Bild\"> <span style=\"color: #000000;\">" + build + "</span>\n"
            data += "      </a><br>\n"

            data += "      <div id=\"file_" + str(f) + "\" class=\"Rand\" style=\"display: none;\">\n"
            data += "<!--- " + str(c) + ". File -->\n"
            if os.path.isdir(p + "/" + computer + "/" + build + "/" + exp_plot_Info): 
              for fi in os.listdir(p + "/" + computer + "/" + build + "/" + exp_plot_Info + "/plots/"):
#	      if fi == exp_file_Info + ".png":
                data += "<a href=\"#\" class=\"link\"><span style=\"color: #FF0000;\">" + fi.strip(".png") + "</span></a><br>\n"
            else:
              data += "<a href=\"#\" class=\"link\"><span style=\"color: #FF0000;\"> No File exist !!! </span></a><br>\n"
            
            data += "<!--- Ende für file -->\n"
            data += "      </div>\n"
           
          data += "    </div>\n"
          
        data += "  </div>\n"
        data += "  </div>\n"
        data += "</div>\n"
              
	data += '''
          </td>
        </tr>
      </table>
    </div>
  </div>
  <div style="float:left; padding:3px; margin:5px  width:200px;"> 
    <div id="arch_ref">
	  '''
        p = "./public_html/archive/ref_" + exp_plot_Info
        if os.path.isdir(p):
	  ref_date = os.listdir(p)[-1]

          p2 = p + "/" + ref_date + "/buildbot/"
	  ref_rev = os.listdir(p2)[-1]

          p3 = p2 + "/" + ref_rev
	  ref_comp = os.listdir(p3)[-1]

          p4 = p3 + "/" + ref_comp
	  ref_build = os.listdir(p4)[-1]
	  
          data += "<h1>Reference Plot (" + ref_date + ")</h1>"
          data += "<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"600\">\n"
	  data += "  <tr>\n    <td>"
          data += "      <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
          data += "        <tr><td><b>Rev. Nr:</b></td><td>" + ref_rev + "</td>\n"
          data += "            <td><b>Computer:</b></td><td>" + ref_comp +"</td>\n"
          data += "            <td><b>Build. Nr:</b></td><td>" + ref_build + "</td>\n"
          data += "            <td> <input type=\"button\" name=\"Text 1\" value=\"Replace\" ></td>\n"
          data += "</tr>\n      </table>\n"
          data += "    </td></tr>\n    <tr><td>\n"
	  data += "      <img src=\"archive/ref_"+exp_plot_Info+"/"+ref_date+"/buildbot/"+ref_rev+"/"+ref_comp+"/"+ref_build+"/"+exp_plot_Info+"/plots/"+exp_file_Info+"\"/>\n"
	  data += "    </td>"
	  data += "  </tr>"
	  data += "</table>"
        else:
          data += "<h1>Reference Plot (xxxx-xx-xx)</h1>"
          data += "<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"600\">\n"
	  data += "  <tr>\n    <td>"
          data += "      <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
          data += "        <tr><td><b>Rev. Nr:</b></td><td> ????</td>\n"
          data += "            <td><b>Computer:</b></td><td> ???? </td>\n"
          data += "            <td><b>Build. Nr:</b></td><td> ???? </td>\n"
          data += "            <td> <input type=\"button\" name=\"Text 1\" value=\"Replace\" ></td>\n"
          data += "</tr>\n      </table>\n"
          data += "    </td></tr>\n    <tr><td>\n"
#	  data += "      <img src=\"archive/ref_" + exp_plot_Info + "/" + ref_date + "/buildbot/" + ref_rev + "/" + ref_comp + "/" + ref_build + "/" + exp_plot_Info + "/plots/" + exp_file_Info + "\" />\n"
	  data += "    </td>"
	  data += "  </tr>"
	  data += "</table>"
	  
	data += '''
    </div>
    <div id="arch_plots">
      <h1>Plot Area</h1>
	  '''

        p = "./public_html/archive/" + svn_Date + "/buildbot/" + svn_Rev + "/"
        data += "  <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
        
        for comp in Archive_Button_Dict['comp']:
#        for computer in os.listdir(p):
	  if os.path.isdir(p + "/" + comp): 
            for build in os.listdir(p + "/" + comp + "/"):
              if os.path.isdir(p + "/" + comp + "/"+ build + "/" + exp_plot_Info): 
	        data += "<tr>"
	        data += "  <td>"
	        data += "    <table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"400\">\n"
#	        data += "      <tr><td>Date:</td><td>" + svn_Date + "</td></tr>\n"
#	        data += "      <tr><td>Rev. Nr:</td><td>" + svn_Rev + "</td></tr>\n"
	        data += "      <tr><td><b>Computer:</b></td><td>" + comp + "</td><td><b>Build. Nr:</b></td><td>" + build + "</td></tr>\n"
#	        data += "      <tr><td>Exp:</td><td>" + exp_plot_Info + "</td></tr>\n"
	        data += "    </table>\n"
                data += "  </td></tr>\n"
                data += "  <tr><td>\n"        
	        data += "    <img src=\"archive/" + svn_Date + "/buildbot/"+ svn_Rev + "/" + comp +"/"+ build + "/" + exp_plot_Info + "/plots/" + exp_file_Info + "\" />\n"
	        data += "  </td>\n"
	        data += "</tr>\n"
        data += "</table>\n"

	data += '''
    </div>
  </div>
        '''

        return data
    
    def get_info(self, request):
        global exp_plot_Info
        global l_nightly
	
        exp_plot_Info="NotSet"
        l_nightly="NotSet"
        
        if "exp" in request.args:
          try:
            exp_plot_Info = request.args["exp"][0]
          except ValueError:
            pass

        if "modus" in request.args:
          try:
            l_nightly = request.args["modus"][0]
          except ValueError:
            pass
	  
        return None

#================================== NEW ================================
