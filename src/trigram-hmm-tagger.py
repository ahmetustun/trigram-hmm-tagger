'''
Created on Mar 10, 2013

@author: arenduchintala
'''
import math
import sys, os, re
from pprint import pprint


unigramCounts = {}
bigramCounts = {}
trigramCounts = {}
emissionCounts = {}
seenObservations = {}
emissions = {}
transitions = {}
pi = {}
s2i = {'I-GENE':2, 'O':3, '*':1, 'STOP':0}
i2s = {}

def populateStateMap(_state):
    return s2i[_state]
    # s2i[_state] = _state
    # return s2i[_state]
    # if (not s2i.has_key(_state)):
    #   s2i[_state] =  len(s2i)
    # print _state, s2i[_state]
    # return s2i[_state]
                
def MLEstimates(countFilePath):
    populateStateMap('*')
    populateStateMap('STOP')
    
    lines = open(countFilePath, 'r').readlines()
    for line in lines:
        line = line.rstrip()
        # print line
        items = line.strip().split(" ")
        if (items[1].strip() == 'WORDTAG'):
            # 8 WORDTAG O localizes
            emissionCounts[(populateStateMap(items[2]) , items[3])] = int(items[0])
        elif (items[1].strip() == '1-GRAM'):
            # 749 1-GRAM I-GENE
            
            unigramCounts[populateStateMap(items[2])] = int(items[0])
            
                
        elif (items[1].strip() == '2-GRAM'):
            # 749 2-GRAM * I-GENE
            bigramCounts[(populateStateMap(items[2]) , populateStateMap(items[3]))] = int(items[0])
        elif (items[1].strip() == '3-GRAM'):
            # 1813 3-GRAM I-GENE O STOP
            trigramCounts[(populateStateMap(items[2]) , populateStateMap(items[3]), populateStateMap(items[4]))] = int(items[0])

    # estimate emission probabilities
    for k, v in emissionCounts.items():
        state = k[0]
        obs = k[1]
        key = (obs , state)
        seenObservations[obs] = 1
        value = math.log(v / float(unigramCounts[state]))
        # emissions[key] = value
        emissions[(obs, state)] = value
    
    # estimate transition probabilities
    for k, v in trigramCounts.items():
        # k = D,N,V 
        state_old = (k[0] , k[1])
        # count of D,N,V over count of D,N
        value = math.log(v / float(bigramCounts[state_old]))
        # probability_key = V/D,N
        probability_key = (k[2] , k[0] , k[1])
        # transitions[probability_key] = value
        transitions[probability_key] = value
    

def getWordPossibleTags(word):    
    states = []
    for state in unigramCounts.keys():
        if emissions.has_key((word , state)):
            states.append(state)
    return states

def oldRareWords(token):
    if len(re.findall(r'[0-9]', token)) > 0:
        return '_RARE_NUMERIC_'
    elif  len(re.findall(r'[A-Z]', token)) == len(token):
            return '_RARE_ALL_CAPS_'
    elif len(re.findall(r'[A-Z]', token[len(token) - 1])) == 1:
            return '_RARE_LAST_CAP_'
    else: 
        return '_RARE_'
        
def replaceRareWord(token):
    return oldRareWords(token)
    # try:
    #    float(token)
    #    return '_RARE_NUMERIC_'
    # except:
    #    
    #    if (token.upper() == token):
    #        return '_RARE_ALL_CAPS_'
    #    elif (token[len(token)-1].upper() == token[len(token)-1]):
    #        return '_RARE_LAST_CAP_'
    #    elif (token[0].upper() == token[0]):
    #        return '_RARE_FIRST_CAP_'
    #    elif (token.endswith('s')):
    #        return '_RARE_END_S_'
    #    else:
    #        return '_RARE_'

    
        
def getViterbiProbability(possible_states_per_word, words):
    pi = {}
    arg_pi = {}
    pi[(0, s2i['*'], s2i['*'])] = 0  # 0,*,*
    arg_pi[(0, s2i['*'], s2i['*'])] = []
    
    for k in range(1, max(possible_states_per_word.keys())+1):
        for v in possible_states_per_word[k]:
            for u in possible_states_per_word[k - 1]:
                prob2bt = {}
                for w in possible_states_per_word[k - 2]:
                    print 'k=', k, 'w=', w, 'u=', u, 'v=', v, words[k]
                    # if (transitions.has_key(v + '|' + w + ',' + u)):
                    if transitions.has_key((v, w , u)):
                        q = transitions[(v , w , u)]
                        print str(v) , '|' , str(w) , ',' + str(v) , '=' , str(q)
                    else:
                        q = float("-inf")                     
                    
                    if emissions.has_key((words[k], v)):
                        e = emissions[(words[k] , v)]
                        print words[k], '|', v , '=' , str(e)
                    else:
                        e = float("-inf")
                   
                    pi_key = str(k - 1) + ',' + str(w) + ',' + str(u)
                    print 'searching pi_key', pi_key
                    
                    p = pi[(k - 1 , w , u)] + q + e
                    bt = list(arg_pi[(k - 1, w , u)])              
                    bt.append(w)
                    prob2bt[p] = bt
                max_p = max(prob2bt.keys())
                print 'max_p =', max_p, 'when k,v,w,u' , k, v, w, u
                pprint(prob2bt)
                max_bt = prob2bt[max_p]
                new_pi_key = (k, u , v)
                pi[new_pi_key] = max_p
                arg_pi[new_pi_key] = max_bt

    # finding max u,v
    k =  max(possible_states_per_word.keys())
    prob2bt = {}
    for v in possible_states_per_word[k]:
        for u in possible_states_per_word[k - 1]:
                print 'last k,u,v,word=', k, u, v ,words[k-1], words[k] 
                if transitions.has_key((s2i['STOP'], u , v)):
                    q = transitions[(s2i['STOP'], u , v)]
                    print s2i['STOP'] , '|' , str(u) , ',' + str(v) , '=' , str(q)
                else:
                    q = float("-inf") 
                p = pi[(k  , u , v)] + q
                bt = list(arg_pi[(k, u , v)])             
                bt.append(u)
                bt.append(v)
                prob2bt[p] = bt;

    max_bt = prob2bt[max(prob2bt.keys())]
    max_p = max(prob2bt.keys())
    print max_p, max_bt
    max_bt.pop(0)
    max_bt.pop(0)
    max_bt = map(lambda bt: i2s[bt], max_bt)
    print "back trace:", max_bt
    print "probability", max_p
    return max_bt



def getPossibleSates(sentence):
    words = {}
    possible_states_per_word = {}
    possible_states_per_word[-1] = [s2i['*']]
    possible_states_per_word[0] = [s2i['*']]
    i = 1
    for word in sentence.split("\n"):
        word = word.strip()
        if (word == ''):
            print 'skipping empty word'
        else:
            if (seenObservations.has_key(word)):
                words[i] = word
                possible_states_per_word[i] = getWordPossibleTags(word)
            else:
                rare_word = replaceRareWord(word)
                words[i] = rare_word
                possible_states_per_word[i] = getWordPossibleTags(rare_word)
            # print word, possible_states_per_word[i]
            i += 1
    
    #possible_states_per_word[i] = [s2i['STOP']]
    return (possible_states_per_word, words)


    
if __name__ == "__main__":

    countFile = '../outputs/gene-rare-cat.count'
    testFile = '../data/gene.dev'
    outputFile = '../outputs/gene-dev-with-state-map.p3.out'
    writer = open(outputFile, 'w')
    MLEstimates(countFile)
    i2s = {v:k for k, v in s2i.items()}  # create a number to state mapping for displaying the state sequence in the end
    print s2i
    print i2s
    testSentences = open(testFile, 'r').read().split("\n\n")
    for i, a_sentence in enumerate(testSentences):
        (possible_states_per_word, words) = getPossibleSates(a_sentence.rstrip())
        # print a_sentence.strip().split("\n")
        if(len(words) > 0):
            tags = getViterbiProbability(possible_states_per_word, words)
        else:
            tags  =[]
        original_words = a_sentence.strip().split('\n')
        for word, tag in zip(original_words, tags):
            # print word, tag
            writer.write(word + ' ' + tag + '\n')
        if(i < len(testSentences) - 1):
            writer.write('\n')    
    writer.flush()
    writer.close()
    os.system("python ../src/eval_gene_tagger.py ../data/gene.key " + outputFile)
    
