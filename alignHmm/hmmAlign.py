import math
import fasta
from context import *
import numpy as np
from hmmmodel 
import numpy as np

class HmmAlign:
     def __init__(self, fastaName, hmmmodel):
        self.seq = fasta.load(fastaName)[1]
        self.hmmmodel = hmmmodel 
        self.states = ['M', 'I', 'D']

    def viterbi(self):
        #searching seq
        seq = self.seq
        stateL = [] 
        v_m = 0 
        v_i = 0 
        v_d = 0 

        for i in range(1, len(seq)):
            aa = seq[i]
            e_m = hmmmodel.getEmitProb(e_i, aa) 
            qxi = 666
            beforeTerm_m = math.log(e_m/qxi) 

            prev_v_m = v_m
            prev_v_i = v_i
            prev_v_d = v_d

            if i < len(seq) -1: # middle states
                amm = hmmmodel.getTransitProb(('M', i), ('M', i+1))
                aim = hmmmodel.getTransitProb(('I', i), ('M', i))
                adm = hmmmodel.getTransitProb(('D', i), ('M', i+1))

                ami = hmmmodel.getTransitProb(('M', i), ('I', i))
                aii = hmmmodel.getTransitProb(('I', i), ('I', i))
                adi = hmmmodel.getTransitProb(('D', i), ('I', i))

                amd = hmmmodel.getTransitProb(('M', i), ('D', i+1))
                aid = hmmmodel.getTransitProb(('I', i), ('D', i))
                add = hmmmodel.getTransitProb(('D', i), ('D', i+1))

                v_m = beforeTerm_m  + max( prev_v_m + amm, prev_v_i + aim, prev_v_d + adm)
                v_i =  max( prev_v_m + ami, prev_v_i + aii, prev_v_d + adi)
                v_d =  max( prev_v_m + amd, prev_v_i + aid, prev_v_d + add)

            else: # last state????
                amm = hmmmodel.getTransitProb(('M', i), ('M', i+1))
                aim = hmmmodel.getTransitProb(('I', i), ('M', i))
                adm = hmmmodel.getTransitProb(('D', i), ('M', i+1))

                ami = hmmmodel.getTransitProb(('M', i), ('I', i))
                aii = hmmmodel.getTransitProb(('I', i), ('I', i))
                adi = hmmmodel.getTransitProb(('D', i), ('I', i))

                amd = hmmmodel.getTransitProb(('M', i), ('D', i+1))
                aid = hmmmodel.getTransitProb(('I', i), ('D', i))
                add = hmmmodel.getTransitProb(('D', i), ('D', i+1))

                v_m =  max( prev_v_m + amm, prev_v_i + aim, prev_v_d + adm)
                v_i =  max( prev_v_m + ami, prev_v_i + aii, prev_v_d + adi)
                v_d =  max( prev_v_m + amd, prev_v_i + aid, prev_v_d + add)


            which_best = np.argmax([v_m, v_i, v_d])
            best_score = np.amax([v_m, v_i, v_d])
            best_state = self.states[which_best]
            stateL.append(best_state)
            
        return best_score, stateL
    
    def retrieveAlignment(self, stateL):
        seq = self.seq 
        alignment_seq = []
        for i in range(len(seq)):
            state = stateL[i]
            if state == 'M':
                alignment_seq.append(seq[i])
            else:
                alignment_seq.append('-')








