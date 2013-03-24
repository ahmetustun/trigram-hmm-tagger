'''
Created on Mar 15, 2013

@author: arenduchintala
'''

if __name__ == "__main__":
    preds = open('../conll/pos.out','r').readlines()
    trues =  open('../conll/pos.key','r').readlines()
    correct = 0
    total = 0
    for ps, ts in zip(preds, trues):
        if ps == ts:
            correct += 1
        else:
            print ps.rstrip(), ts.rstrip()
        total +=1
    
    print float(correct)/total
        
    