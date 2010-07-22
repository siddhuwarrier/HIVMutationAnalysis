'''
Created on 22 Jul 2010

@author: siddhu
'''
from google.appengine.ext import db

class FastaFileDb(db.Model):
    author = db.UserProperty()
    #the fasta file itself
    fasta = db.StringProperty(multiline = True)
    
    
    