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
	
      data  = '''
<h1>ICON Buildbot</h1>
Buildbot starts every night a test suite on a set of builders differing in machines, compilers 
and parallelization setup. A standard Buildbot test unpacks the model code and scripts from the  
repositors, configures and compiles the codes for the selected builder, and runs a series of 
standard test experiments and related post-processings. Figures resulting from the 
post-processings are archived and are uploaded to web pages for comparison with reference 
figures.<p>
      '''
	
      data += '''
<h2>Status results</h2>
<table border="1" cellspacing="0" cellpadding="2">
   <thead>
     <tr valign="baseline">
       <td><i>Display</i></td>
       <td><i>What it shows</i></td>
     </tr>
   </thead>
   <tfoot></tfoot>
   <tbody>
     <tr valign="baseline">
       <td><a href="waterfall?show_events=false&show_time=604800">Waterfall</a></td>
       <td>Status of each step of the Buildbot tests for the past 7 days</td>
     </tr>

     <tr valign="baseline">
       <td><a href="grid?width=10">Grid</a></td>
       <td>Status of the complete Builbot test for the 10 latest tested revisions</td>
     </tr>

     <tr valign="baseline">
       <td><a href="one_box_per_builder">Latest build</a></td>
       <td>Status of the complete Builbot test for the latest test on each builder</td>
     </tr>

     <tr valign="baseline">
       <td><a href="one_line_per_build">List</a></td>
       <td>Status of the complete Builbot test for the latest 15 tests</td>
     </tr>

     <tr valign="baseline">
       <td><a href="buildslaves">Build slaves</a></td>
       <td>Information on machines ("slaves"), builders, etc.</td>
     </tr>
   </tbody>
   </table>
      '''

      data += '''
<h2>Eperiment results</h2>
<table border="1" cellspacing="0" cellpadding="2">
  <thead>
    <tr valign="baseline">
      <td width="300"><i>Experiment (exp)</i></td>
      <td width="380"><i>Description</i></td>
      <td width="200"><i>Model</i></td>
      <td width="120"><i>Grid</i></td>
    </tr>
  </thead>
  <tfoot>
  </tfoot>
  <tbody>
      '''
        
# Build HTML Table info

      for exp in ExpList:
        for e in ExpListInfo:
          if exp == e.get('ExpName'):
	    break
	      
 	data += '''
<tr valign="baseline">
        '''
        data += "<td><a href=\"plot?exp=" + e.get('ExpName') + "&modus=nightly\">" + e.get('ExpName') + "</a></td>"
        data += "<td>" + e.get('Description') + "</td>"
        data += "<td>" + e.get('Model') + "</td>"
        data += "<td>" + e.get('Grid') + "</td>"
        data += "</tr>"
            

      data += '''
 	  </tbody>
	  </table>
      '''

      data += '''
          <h3>About BuildBot</h3>
	  <ul>
	    <li><a href="about">Buildbot version used here</a></li>
	    <li>Buildbot home page: <a href="http://trac.buildbot.net" target="_blank">http://trac.buildbot.net</a></li>
	  </ul>
      '''

      return data
  
   
    def get_info(self, request):
        global change_Info
	
        change_Info="NotSet"
#        exp_Info="NotSet"
        
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
	print "WS: ======"
	print data
	print "WS: ======"
        return data
#================================== NEW ================================
