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
  was.src="gif/open.gif";
 }
 else {
  document.getElementById(das).style.display='none';
   was.src='gif/closed.gif';
  }
}
 </script>
        '''
        return h
						
    def body(self, request):
        global exp_plot_Info
        global l_nightly
	self.get_info(request)
	
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
         
	data = '''
  <div style="float:left; padding:3px; margin:5px width:350px;">
    <div id="arch_menu" >
      <h1>ICON Buildbot Archive</h1>
      <table border=""1 cellspacing="0" cellpadding="0" align="left" width="400">
        <colgroup>
          <col width="100">
          <col width="250">
        </colgroup>

        <tr>
          <td><b>Date</b> <br> from/to</td>
          <td>
	    <form name=date_form>
	    
	'''

	data += "<select  disabled=\"disabled\">\n"
        data += "<option selected>" + tmp_Date.split('-')[0] + "</option></select>\n"
        data += "</select>\n"
	
	data += "<select  disabled=\"disabled\">\n"
        data += "<option selected>" + tmp_Date.split('-')[1] + "</option></select>"
        data += "</select>\n"
	
	data += "<select  disabled=\"disabled\">\n"
        data += "<option selected>" + tmp_Date.split('-')[2].split(' ')[0] + "</option></select>\n"
        data += "</select>\n"
	
	data += '''
	    </form>
	    <form name=date_form>
	'''

	data += "<select  disabled=\"disabled\">\n"
        data += "<option selected>" + tmp_Date.split('-')[0] + "</option></select>"
        data += "</select>\n"
	
	data += "<select  disabled=\"disabled\">\n"
        data += "<option selected>" + tmp_Date.split('-')[1] + "</option></select>"
        data += "</select>\n"
	
	data += "<select  disabled=\"disabled\">\n"
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
	
	exp_file_Info = "test_hat_jww_R2B04L31_DIV850"
	data += '''
	      </select>
	    </form>
          </td>
        </tr>
        <tr>
          <td><b>File</b></td>
            <td>
	      <form name=file_form>
		<select name=mytextarea disabled="disabled">
	      '''
	data += "<option selected>" + exp_file_Info + "</option>"
	data += '''
	        </select>
	      </form>
            </td>
          </td>
        </tr>
	'''
        
        data += "<tr><td><b>Search</b><td>\n"
        data += "  <input type=\"button\" name=\"Text 1\" value=\"OK\" disabled=\"disabled\">\n"
        data += "</td></tr>\n"
        
	data += '''
      </table>
    </div>
    <div id="arch_file_name">
      <h1>Available Plot List</h1>
      <table border="0" cellspacing="0" cellpadding="0" align="left" width="200">
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
            for fi in os.listdir(p + "/" + computer + "/" + build + "/" + exp_plot_Info + "/plots/"):
#	      if fi == exp_file_Info + ".png":
              data += "<a href=\"#\" class=\"link\"><span style=\"color: #FF0000;\">" + fi + "</span></a><br>\n"
            data += "<!--- Ende für file -->\n"
            data += "      </div>\n"
           
          data += "<!--- Ende für Builds build -->\n"
          data += "    </div>\n"
          
        data += "<!--- Ende für Computer -->\n"
        data += "  </div>\n"
        data += "<!--- Ende für Date -->\n"
        data += "  </div>\n"
        data += "<!--- Ende für Rev -->\n"
        data += "</div>\n"
              
	data += '''
          </td>
        </tr>
      </table>
    </div>
  </div>

  <div style="float:left; padding:3px; margin:5px  width:200px;"> 
    <div id="arch_ref">
      <h1>Reference Plot</h1>
      <table border="0" cellspacing="0" cellpadding="0" align="left" width="600">
	<tr>
	  <td>
            <table border="0" cellspacing="0" cellpadding="0" align="left" width="200">
              <tr>   <td>Date:</td>      <td>2011-06-20</td>  </tr>
              <tr>   <td>Rev. Nr:</td>   <td>4725</td>        </tr>
              <tr>   <td>Computer:</td>  <td>TORNADO_pgi</td>        </tr>
              <tr>   <td>Build. Nr:</td>   <td>582</td>        </tr>
              <tr>   <td>Exp:</td>       <td>test_hat_jww</td>        </tr>
            </table>
          </td>
          <td>
	  '''
        
	data += "<img src=\"archive/" + svn_Date + "/buildbot/"+ svn_Rev +"/TORNADO_sun/1160/" + exp_plot_Info + "/plots/" + exp_file_Info + ".png\" />\n"
	data += '''
          </td>
        </tr>
      </table>
    
    </div>
    <div id="arch_plots">
      <h1>Plot Area</h1>
      <table border="0" cellspacing="0" cellpadding="0" align="left" width="400">
        <!--- Plot line one --->
          <tr>
	    <td>
	      <table border="0" cellspacing="0" cellpadding="0" align="left" width="200">
	        <tr>   <td>Date:</td>      <td>2011-06-20</td>  </tr>
		<tr>   <td>Rev. Nr:</td>   <td>4725</td>        </tr>
		<tr>   <td>Computer:</td>  <td>TORNADO_pgi</td>        </tr>
		<tr>   <td>Build. Nr:</td>   <td>582</td>        </tr>
		<tr>   <td>Exp:</td>       <td>test_hat_jww</td>        </tr>
	      </table>
            </td>
            <td>
	  '''
        
	data += "<img src=\"a/" + svn_Date + "/buildbot/"+ svn_Rev +"/TORNADO_sun/1160/" + exp_plot_Info + "/plots/" + exp_file_Info + ".png\" />\n"
	data += '''
	    </td>
	  </tr>


	  <!--- Plot line 2 --->
	  <tr>
	    <td>
              <table border="0" cellspacing="0" cellpadding="0" align="left" width="200">
		<tr>   <td>Date:</td>      <td>2011-06-20</td>  </tr>
	        <tr>   <td>Rev. Nr:</td>   <td>4725</td>        </tr>
                <tr>   <td>Computer:</td>  <td>BLIZZ_yMyO</td>        </tr>
                <tr>   <td>Build. Nr:</td>   <td>557</td>        </tr>
                <tr>   <td>Exp:</td>       <td>test_hat_jww</td>        </tr>
              </table>
           </td>
           <td>
	      <img src="archive/2012-04-04/buildbot/8609/TORNADO_sun/1068/test_hat_jww/plots/test_hat_jww_R2B04L31_DIV850_day09.png" />
	   </td>
	 </tr>

	    <!--- Plot line 3 --->

	 <tr>
	   <td>
	     <table border="0" cellspacing="0" cellpadding="0" align="left" width="200">
	       <tr>   <td>Date:</td>      <td>2011-06-20</td>  </tr>
	       <tr>   <td>Rev. Nr:</td>   <td>4725</td>        </tr>
               <tr>   <td>Computer:</td>  <td>BLIZZ_nMnO</td>        </tr>
               <tr>   <td>Build. Nr:</td>   <td>566</td>        </tr>
               <tr>   <td>Exp:</td>       <td>test_hat_jww</td>        </tr>
             </table>
           </td>
           <td>
	      <img src="archive/2012-04-04/buildbot/8609/TORNADO_sun/1068/test_hat_jww/plots/test_hat_jww_R2B04L31_DIV850_day09.png" />
	   </td>
	 </tr>

	    <!--- Plot line 4 --->

	 <tr>
	   <td>
	     <table border="0" cellspacing="0" cellpadding="0" align="left" width="200">
	       <tr>   <td>Date:</td>      <td>2011-06-20</td>  </tr>
	       <tr>   <td>Rev. Nr:</td>   <td>4725</td>        </tr>
               <tr>  <td>Computer:</td>  <td>TORNADO_nag</td>        </tr>
               <tr>   <td>Build. Nr:</td>   <td>605</td>        </tr>
               <tr>   <td>Exp:</td>       <td>test_hat_jww</td>        </tr>
             </table>
          </td>
          <td>
	      <img src="archive/2012-04-04/buildbot/8609/TORNADO_sun/1068/test_hat_jww/plots/test_hat_jww_R2B04L31_DIV850_day09.png" />
	  </td>
	</tr>
	    <!--- Plot line 5 --->
	<tr>
	  <td>
	    <table border="0" cellspacing="0" cellpadding="0" align="left" width="200">
	      <tr>   <td>Date:</td>      <td>2011-06-20</td>  </tr>
	      <tr>   <td>Rev. Nr:</td>   <td>4725</td>        </tr>
              <tr>   <td>Computer:</td>  <td>BLIZZ_yMnO</td>        </tr>
              <tr>   <td>Build. Nr:</td>   <td>731</td>        </tr>
              <tr>   <td>Exp:</td>       <td>test_hat_jww</td>        </tr>
            </table>
          </td>
          <td>
	      <img src="archive/2012-04-04/buildbot/8609/TORNADO_sun/1068/test_hat_jww/plots/test_hat_jww_R2B04L31_DIV850_day09.png" />
	  </td>
	</tr>

	    <!--- Plot line 6 --->
	<tr>
	  <td>
	    <table border="0" cellspacing="0" cellpadding="0" align="left" width="200">
	      <tr>   <td>Date:</td>      <td>2011-06-20</td>  </tr>
	      <tr>   <td>Rev. Nr:</td>   <td>4725</td>        </tr>
	      <tr>   <td>Computer:</td>  <td>BLIZZ_nMyO</td>        </tr>
              <tr>   <td>Build. Nr:</td>   <td>558</td>        </tr>
              <tr>   <td>Exp:</td>       <td>test_hat_jww</td>        </tr>
            </table>
          </td>
          <td>
              <img src="archive/2012-04-04/buildbot/8609/TORNADO_sun/1068/test_hat_jww/plots/test_hat_jww_R2B04L31_DIV850_day09.png" />
          </td>
	</tr>        
    </table>

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
