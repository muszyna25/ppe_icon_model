#================================== NEW ================================

from buildbot.status.web.base import HtmlResource
from twisted.web.util import Redirect
import os, time
from shutil import copytree
from buildbot import util

class MainPage(HtmlResource):
    title = "Replace Reference Buildbot"

    def head(self, request):
        h='''
<link href="archive.css" rel="stylesheet" type="text/css" />
        '''
        return h

    def body(self, req):
      
#=========================================================

#=========================================================

      def check_rev_exist(exp,rev):
        st = False
        c_list = []
        d_list = []
        p = "public_html/archive/"
        for d in os.listdir(p):
	  if os.path.isdir(p + d):
	    p_d = p + d + "/buildbot/"
            for r in os.listdir(p_d):
              if r == rev:
	        p_r = p_d + r + "/"
                for c in os.listdir(p_r):
	          p_c = p_r + c + "/"
                  for b in os.listdir(p_c):
	            p_b = p_c + b + "/"
	            if os.path.isdir(p_b + exp):
		      st = True
                      d_list.append(d)
                      c_list.append(c)
                        
        d_list = list(sorted(d_list))                
        c_list = list(sorted(set(c_list)))                
        return (st,d_list,c_list)
#=========================================================

#=========================================================
      def create_BuildList(exp,rev,d_list,comp):
        st = False
        b_list = []
        p = "public_html/archive/"
        for d in d_list:
	  p_d = p + d + "/buildbot/"
          p_r = p_d + rev + "/"
	  if os.path.isdir(p_r + comp):
	    p_c = p_r + comp + "/"
            for b in os.listdir(p_c):
	      p_b = p_c + b + "/"
	      if os.path.isdir(p_b + exp):
	        st = True
                b_list.append(b)
                        
        b_list = list(sorted(set(b_list)))                
        return (st,b_list)

#=========================================================

#=========================================================
      def get_Date(exp,rev,comp,build):
        st = False
        b_list = []
        pr = "public_html/archive/"
        for d in os.listdir(pr):
#        for d in d_list:
#         public_html/archive/yyyy-mm-dd/buildbot/
	  p = pr + d + "/buildbot/"
#         public_html/archive/yyyy-mm-dd/buildbot/rrrr/
	  p += rev + "/"
#         public_html/archive/yyyy-mm-dd/buildbot/rrrr/ccccc
	  p += comp + "/"
#         public_html/archive/yyyy-mm-dd/buildbot/rrrr/ccccc/bbbb
	  p += build + "/"
#         public_html/archive/yyyy-mm-dd/buildbot/rrrr/ccccc/bbbb/eeee
	  p += exp + "/"
	  if os.path.isdir(p):
	    return d
	    
	return "yyyy-mm-dd"

#-------------------------------------------------------------
      e  = "Experement Name"
      le = False
      lr = False
      lc = False
      lb = False
      #ls = False

      l_use_build = False
      l_use_comp  = False
      
      if "exp" in req.args:
        try:
          EXP = req.args["exp"][0]
          le = True
        except ValueError:
          pass

      if "rev" in req.args:
        try:
          REV = req.args["rev"][0]
          lr = True
        except ValueError:
          pass

      if "comp" in req.args:
        try:
          COMP = req.args["comp"][0]
          lc = True
        except ValueError:
          pass

      if "build" in req.args:
        try:
          BUILD = req.args["build"][0]
          lb = True
        except ValueError:
          pass
	
      #if "button" in req.args:
        #try:
          #CONTROL = req.args["button"][0]
          #ls = True
        #except ValueError:
          #pass
	
      #if ls:
        #if CONTROL == "cancel":
	  #return Redirect("../plot")
	
        #if CONTROL == "ok":
	  #return Redirect("../plot")
	
      print "===== WS ===="
      print req
      print req.args
      print "===== WS ===="

      data = "<div id=\"div_ref\"  style=\"float:left; padding:3px; margin:5px width:500px;\">"
      data += "<h1>Replace Reference Plot</h1>\n"
      
      if not le:
        data += "<br><br><br><h2>Pleace use ?exp=exp_name \n"
        return data
	
      data += "<h2>" + EXP + "</h1>\n"
      data += "<table border=\"0\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\" width=\"300\">"

# Set Revision Info

      data += "  <tr>"
      data += "    <td>Revision:</td>"
      if not lr:
        data += "  <td>"
        data += "    <form name=\"replace_save\" method=\"POST\" action=\"/reference?exp=" + EXP + "\"" + " class=\"command replace\">\n"
        data += "      <input name=\"rev\" type=\"text\" size=\"5\" maxlength=\"10\">\n"
        data += "      <input type=\"submit\" name=\"Revision_button\" value=\">\" >\n"
        data += "    <form>"
        data += "  </td>"
      else:
	data += "<td>" + REV + "</td>"
        l_use_comp  = True
      data += "  </tr>"
      
      c_list = ["----"]
      b_list = ["----"]
      
#      if l_use_comp and not lc:
      if l_use_comp:
	st,d_list,c_list = check_rev_exist(EXP,REV)
	if not st:
          data = "<div id=\"div_ref\" >"
          data += "<h1>Replace Reference Plot</h1>\n"
          data += "Revision " + REV + " contains no runs for experiment "+ EXP 
          data += "</div>"
          return data

	

# Set Computer Info

      data += "  <tr>\n"
      data += "    <td>Computer:</td>\n"
      if not lc:
        data += "  <td>\n"
        data += "    <form name=\"reference_comp\" method=\"POST\" action=\"/reference?exp=" + EXP
        if not l_use_comp:
          data += "\" class=\"command replace\">\n"
          data += "      <select name=\"comp\" disabled=\"disabled\">\n"
        else:
          data += "&rev=" + REV
          data += "\" class=\"command replace\">\n"
          data += "      <select name=\"comp\">\n"
        for c in c_list:
          data += "        <option>" + c + "</option>\n"
        data += "      </select>"
        if l_use_comp:
          data += "      <input type=\"submit\" name=\"Computer_button\" value=\">\" >\n"
        data += "    <form>\n"
        data += "  </td>"
      else:
	data += "    <td>" + COMP + "</td>"
        l_use_build = True
      data += "    </tr>"

# Set Build Info

      if l_use_build:
	st,b_list = create_BuildList(EXP,REV,d_list,COMP)
	if not st:
          data = "<div id=\"div_ref\" >"
          data += "<h1>Replace Reference Plot</h1>\n"
          data += "Revision " + REV + " contains no runs for experiment "+ EXP 
          data += "</div>"
          return data

      data += "  <tr>\n"
      data += "    <td>Build Nr.:</td>\n"
      if not lb:
        data += "  <td>"
        data += "    <form name=\"reference_build\" method=\"POST\" action=\"/reference?exp=" + EXP
        if not l_use_build:
#          data += "&rev=" + REV 
          data += "\" class=\"command replace\">\n"
          data += "      <select name=\"build\" disabled=\"disabled\">\n"
        else:
          data += "&rev=" + REV
          data += "&comp=" + COMP
          data += "\" class=\"command replace\">\n"
          data += "      <select name=\"build\">\n"
        for b in b_list:
          data += "        <option  selected>" + b + "</option>\n"
        if l_use_build:
          data += "      <input type=\"submit\" name=\"Computer_button\" value=\">\" >\n"
        data += "      </select>\n"
        data += "    <form>\n"
        data += "  </td>"
      else:
	data += "<td>" + BUILD + "</td>"

      data += "  </tr>"

      if le and lr and lc and lb:
	DATE = get_Date(EXP,REV,COMP,BUILD)
        data += "  <tr>"
        data += "    <form name=\"replace_save\" method=\"POST\" action=\"/reference/save_cancel\"  class=\"command replace\">\n"
        
        data += "      <input type=\"hidden\" name=\"date\" value=\""+ DATE + "\">\n"
        data += "      <input type=\"hidden\" name=\"exp\" value=\""+ EXP + "\">\n"
        data += "      <input type=\"hidden\" name=\"rev\" value=\""+ REV + "\">\n"
        data += "      <input type=\"hidden\" name=\"comp\" value=\""+ COMP + "\">\n"
        data += "      <input type=\"hidden\" name=\"build\" value=\""+ BUILD + "\">\n"
        data += "      <td><input type=\"submit\" name=\"button\" value=\"ok\" ></td>\n"
        data += "      <td><input type=\"submit\" name=\"button\" value=\"cancel\" ></td>\n"
        data += "    <form>"
        data += "  </tr>"
      data += "  </table>"
      
      data += "</div>"	
      if le and lr and lc and lb:
        data += "<div>"	
	DATE = get_Date(EXP,REV,COMP,BUILD)
	p = "public_html/archive/" + DATE + "/buildbot/" + REV + "/" + COMP + "/" + BUILD + "/" + EXP + "/plots/"
	for f in os.listdir(p):
          data += "<br>\n"
	  data += "<img src=\"archive/"+DATE+"/buildbot/"+REV+"/"+COMP+"/"+BUILD+"/"+EXP+"/plots/" + f + "\"/>\n"
	  
        data += "</div>"	
      return data


    def save_cancel(self, req):
      print "===== WS ===="
      print req
      print req.args
      print "===== WS ===="
      
      if "button" in req.args:
        try:
          CONTROL = req.args["button"][0]
        except ValueError:
          pass
	
      if "exp" in req.args:
        try:
          EXP = req.args["exp"][0]
        except ValueError:
          pass
	
      if "date" in req.args:
        try:
          DATE = req.args["date"][0]
        except ValueError:
          pass
	
      if "rev" in req.args:
        try:
          REV = req.args["rev"][0]
        except ValueError:
          pass
	
      if "comp" in req.args:
        try:
          COMP = req.args["comp"][0]
        except ValueError:
          pass
	
      if "build" in req.args:
        try:
          BUILD = req.args["build"][0]
        except ValueError:
          pass
	
      if CONTROL == "cancel":
	l = "../plot?exp=" + EXP + "&modus=nightly&status=cancel" 
        return Redirect(l)
	
      if CONTROL == "ok":
	save_time = time.strftime("%Y-%m-%d-%H-%M", time.localtime(util.now()))
	tmpname   = "/" + DATE + "/buildbot"
	tmpname  += "/" + REV
	tmpname  += "/" + COMP
	tmpname  += "/" + BUILD
	
	srcname   = "public_html/archive" + tmpname + "/" + EXP
	
	dstname   = "public_html/reference_plot/"
	dstname  += EXP + "/" + save_time + tmpname
#	dstname  += EXP + "/2012-06-22-19-05/" + tmpname
#	l_dstname = dstname.split('/') 
        print "===== OK ===="
        print srcname
        print dstname
#        print l_dstname
        
        if not os.path.exists(dstname):
          print "Create:" + dstname 
          os.makedirs(dstname)
        else:
          print "OK:" + dstname 

        print "===== OK ===="

	copytree(srcname,dstname + "/"  + EXP)
	l = "../plot?exp=" + EXP + "&modus=nightly&status=ok" 
	return Redirect(l)
	
      return Redirect("../plot")
			      
    def default(self, req):
        return Redirect("home")
							
    def getChild(self, path, req):
        if path == "save_cancel":
	  print "WS: save_cancel"
	  return self.save_cancel(req)
#	else:
#	  print "WS: else"
#	  return self.default(req)

# ================================== NEW ================================
 