'''
Created on 22 Jul 2010

@author: siddhu
'''
import logging
import os.path
from firepython.middleware import FirePythonWSGI

from google.appengine.ext.webapp.util import run_wsgi_app
from google.appengine.ext import webapp
from google.appengine.api import users
from google.appengine.ext.webapp import template

from CodonAnalysis.PositionalAnalysis import FastaAnalysis

###Class for the welcome page which should
# be visible only if user is logged in.
class WelcomePage(webapp.RequestHandler):
    #set up logging
    def get(self):
        logging.getLogger().setLevel(logging.DEBUG)
        logging.debug("Welcome page: get()")
        
        template_args = {
                         'user': users.get_current_user()
                         }
        path = os.path.join(os.path.dirname(__file__), 'templates', 'index.html')
        
        self.response.out.write(template.render(path, template_args))

class ProcessProtSeqPage(webapp.RequestHandler):
    def post(self):
        if self.request.POST.get('fastaFile') == None:
            logging.debug("Redirecting...")
            
        logging.getLogger().setLevel(logging.DEBUG)
        
        #convert the list of integers specified into a pythonic list
        posns = [int(i) for i in self.request.POST.get('codonPositions').split(",") if i.isdigit()]
        
        successfulExecution = True
        selectedProteins = None
        try:
            #perform the actual codon position analysis
            fastaAnalysis = FastaAnalysis(self.request.POST.get('fastaFile').file, posns)
            selectedProteins = fastaAnalysis.createProteinSequences() #selected proteins
        except Exception, error:
            logging.error("Exception occured: " + error.__str__())
            successfulExecution = False
            
        path = os.path.join(os.path.dirname(__file__), 'templates', 'processFasta.html')
        
        templateArgs = {
                        'positions': posns,
                        'selectedProteins': selectedProteins,
                        'success': successfulExecution
                        }
        
        self.response.out.write(template.render(path, templateArgs))
        #seqIterator = SeqIO.parse(self.request.POST.get('myfile').file, "fasta")

        
def main():
    #create the WSGI application object to redirect to the appropriate request
    #handlers 
    application = webapp.WSGIApplication(
                                     [
                                      ('/', WelcomePage),
                                      ('/process', ProcessProtSeqPage)
                                      ],
                                     debug=True)
    #start the webapp
    run_wsgi_app(FirePythonWSGI(application))
    
if __name__ == "__main__":
    main()
