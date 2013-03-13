'''
Created on Mar 11, 2013

@author: arenduchintala
'''

class ViterbiState(object):
    '''
    classdocs
    '''


    def __init__(self, k, maxargs,p):
        '''
        Constructor
        '''
        self.k = k
        self.maxargs = maxargs
        self.p = p
    
    def disp(self):
        print self.k, self.p, self.maxargs
        