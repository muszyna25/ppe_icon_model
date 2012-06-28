#================================== NEW ================================

from buildbot.status.web.base import HtmlResource
#import mainwork
import os

class HomePage(HtmlResource):
    title = "ICON Buildbot"
#    global test_Info=""
#    global change_Info=""
#    global exp_Info=""

    def head(self, request):
        h='''
<link href="archive.css" rel="stylesheet" type="text/css" />
        '''
        return h

    def main_page(self, request):
      global change_Info
      
# Create a emty Experiment List type

      ExpList = []
      ExpListInfo = []

# Read all files in the "icon_run" directory and write the filename which starts with "exp.test_" into the
# Experiment List type

      for filename in os.listdir("icon_run"):
	if filename.find("exp.test_") == 0:
          ExpDict = {"Description" : " ","Model" : " ","Grid" : " "}
          ExpDict["ExpName"] = filename.lstrip("exp.")
            
          file = open("icon_run/" + filename)
          while 1:
            line = file.readline()
            if not line:
              break
             
            if line.find("_Description_") >= 0:
	      tmp,Info = line.split("_bb_table_Description_",1)
              Info.replace("\n","")
              ExpDict["Description"] = Info.strip(" ")
                
            if line.find("_Model_") >= 0:
	      tmp,Info = line.split("_bb_table_Model_",1)
              Info.replace("\n","")
              ExpDict["Model"] = Info.strip(" ")
                
            if line.find("_Grid_") >= 0:
	      tmp,Info = line.split("_bb_table_Grid_",1)
              Info.replace("\n","")
              ExpDict["Grid"] = Info.strip(" ")
              
          file.close
          ExpList.append(filename.lstrip("exp."))
          ExpListInfo.append(ExpDict)
            
            
# Sort the Experiment List type

      ExpList.sort()
	
# Write the upper part of the main page

      data  = '''
<h1>ICON Buildbot</h1>
Buildbot starts every night a test suite on a set of builders differing in machines, 
compilers and parallelization setup. A standard Buildbot test unpacks the latest model 
from the trunk/icon-dev branch of the repository, configures and compiles the codes, 
and runs a series of test experiments and related post-processings. Figures resulting 
from the post-processing are archived and uploaded to web pages for comparison with 
reference figures.<p>
      '''
	
# Write the Status result part to the main page

      data += '''
<h2>Buildbot displays</h2>
<table border="1" cellspacing="5px" cellpadding="2" >
   <thead>
     <tr valign="baseline">
       <td><i>Display</i></td>
       <td><i>Information</i></td>
     </tr>
   </thead>
   <tfoot></tfoot>
   <tbody>
     <tr valign="baseline">
       <td style=\"text-align:left\"><a href="waterfall>Waterfall</a></td>
       <td style=\"text-align:left\">For each step of the test suite, for each builder</td>
     </tr>

     <tr valign="baseline">
       <td style=\"text-align:left\"><a href="grid">Grid</a></td>
       <td style=\"text-align:left\">For the full test suite, for each builder and several revisions</td>
     </tr>

     <tr valign="baseline">
       <td style=\"text-align:left\"><a href="one_box_per_builder">Latest build</a></td>
       <td style=\"text-align:left\">For the full test suite, the latest test on each builder</td>
     </tr>

     <tr valign="baseline">
       <td style=\"text-align:left\"><a href="one_line_per_build">List</a></td>
       <td style=\"text-align:left\">For the full testsuite, the latest 15 completed tests</td>
     </tr>

     <tr valign="baseline">
       <td style=\"text-align:left\"><a href="buildslaves">Build slaves</a></td>
       <td style=\"text-align:left\">Machines and builders</td>
     </tr>
   </tbody>
   </table>
      '''

# Create the experiment List as a table

      
      data += "<h2>Latest nightly test results</h2>\n"
      data += "<table border=\"1\" cellspacing=\"5px\" cellpadding=\"2\">\n"
      data += "  <thead>\n"
      data += "    <tr valign=\"baseline\">\n"
      data += "      <td width=\"300\"><i>Experiment (exp)</i></td>\n"
      data += "      <td width=\"380\"><i>Description</i></td>\n"
      data += "      <td width=\"200\"><i>Model</i></td>\n"
      data += "      <td width=\"120\"><i>Grid</i></td>\n"
      data += "    </tr>\n  </thead>"
      data += "  <tfoot></tfoot>\n"

         
# Build HTML Table info

      data += "  <tbody>\n"
      for exp in ExpList:
        for e in ExpListInfo:
          if exp == e.get('ExpName'):
	    break
	      
 	data += "    <tr valign=\"baseline\">\n"
        data += "      <td style=\"text-align:left\"><a href=\"plot?exp=" + e.get('ExpName') + "&modus=nightly\">" + e.get('ExpName') + "</a></td>\n"
        data += "      <td style=\"text-align:left\">" + e.get('Description') + "</td>\n"
        data += "      <td>" + e.get('Model') + "</td>\n"
        data += "      <td>" + e.get('Grid') + "</td>\n"
        data += "   </tr>"
      
      data += "  </tbody>\n</table>\n"

      data += "<h3>About BuildBot</h3>\n"
      data += "<ul>\n"
      data += "  <li><a href=\"about\">Buildbot version used here</a></li>\n"
      data += "  <li>Buildbot home page: <a href=\"http://trac.buildbot.net\" target=\"_blank\">http://trac.buildbot.net</a></li>\n"
      data += "</ul>\n"

      return data
  
   
    def get_info(self, request):
        global change_Info
	
        change_Info="NotSet"
        
	if "status" in request.args:
          try:
            change_Info = request.args["status"][0]
          except ValueError:
            pass
	
        return None
    
    def body(self, request):
        global change_Info
	self.get_info(request)
        data = self.main_page(request)	
        return data
#================================== NEW ================================
