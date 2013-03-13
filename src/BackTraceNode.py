'''
Created on Mar 12, 2013

@author: arenduchintala
'''

class BackTraceNode(object):
    '''
    classdocs
    '''


    def __init__(self,_state, _p):
        '''
        Constructor
        '''
        self.w = _state
        self.parent = None
        self.prob = _p
        self.children = []
    
    def addChild(self, btNode):
        btNode.parent = self
        self.children.append(btNode)
    
    def disp(self):
        print self.w, self.prob
        