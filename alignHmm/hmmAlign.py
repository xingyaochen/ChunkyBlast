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
        states = self.hmmmodel.matchStates
        #matrices have j rows ,i columns where j range from 0 to len(states)+1, i ranges from 0 to len(seq)-1
        v_m, v_i, v_d = [[None]*(len(seq)+1) for i in range(len(states)+2)], [[None]*len(seq) for i in range(len(states)+2)], [[None]*len(seq) for i in range(len(states)+2)]
        #initialize:

        #initialize M matrix
        v_m[0][0]={"prob": 0, "prev": None}
        for i in range(1,len(seq)+1):
            #begin state M0 cannot match to sequence nucleotide positions 
            v_m[0][i]={"prob": float('-inf'), "prev": None}
        for j in range(1,len(states)+2):
            #can't match fake position to any matching state
            v_m[j][0]={"prob": float('-inf'), "prev": None}
        #end initializing M matrix

        #initialize I matrix:
        for j in range(len(states)+2):
            #can't match fake position to any insertion state
            v_i[j][0]={"prob": float('-inf'), "prev": None}
        #end initialize I matrix

        #initialize D matrix:
        for i in range(1,len(seq)+1):
            #D0 doesn't exit
            v_m[0][i]={"prob": float('-inf'), "prev": None}
        for j in range(len(states)+2):
            #can't match fake position to any insertion state
            v_d[j][0]={"prob": float('-inf'), "prev": None}
        #end initializing D matrix
        qxi=1
        for i in range(1, len(seq)+1):
            aa = seq[i]
            for j in range(1, len(states)+2):
                #already logged
                e_m = hmmmodel.getEmitProb(st, aa) 

                amm = hmmmodel.getTransitProb(('M', j-1), ('M', j))
                aim = hmmmodel.getTransitProb(('I', j-1), ('M', j))
                adm = hmmmodel.getTransitProb(('D', j-1), ('M', j))

                ami = hmmmodel.getTransitProb(('M', j), ('I', j))
                aii = hmmmodel.getTransitProb(('I', j), ('I', j))
                adi = hmmmodel.getTransitProb(('D', j), ('I', j))

                amd = hmmmodel.getTransitProb(('M', j-1), ('D', j))
                aid = hmmmodel.getTransitProb(('I', i-1), ('D', j))
                add = hmmmodel.getTransitProb(('D', j-1), ('D', j))

                # for just v_m 

                m_threeProbs = [v_m[j-1][i-1]["prob"] + amm , v_i[j-1][i-1]["prob"] + aim, v_d[j-1][i-1]["prob"] + adm] 
                i_threeProbs = [v_m[j][i-1]["prob"] + ami , v_i[j][i-1]["prob"] + aii, v_d[j][i-1]["prob"] + adi]
                d_threeProbs = [v_m[j-1][i]["prob"] + amd , v_i[j-1][i]["prob"] + aid, v_d[j-1][i]["prob"] + add]

                max_prob_m = np.amax(m_threeProbs)
                max_prob_i = np.amax(i_threeProbs)
                max_prob_d = np.amax(d_threeProbs)
                prev_state_m = (self.states[np.argmax(m_threeProbs)], j-1, i-1)
                prev_state_i = (self.states[np.argmax(i_threeProbs)], j, i-1)
                prev_state_d = (self.states[np.argmax(d_threeProbs)], j-1, i)
                
                v_m[i][j] = {"prob": math.log(e_m / qxi)  + max_prob_m, "prev": prev_state_m}
                v_i[i][j] = {"prob": max_prob_i, "prev": prev_state_i}
                v_d[i][j] = {"prob": max_prob_d, "prev": prev_state_d}


           
            
        bestScore=v_m[len(states)+1][len(seq)]
        return bestScore
            
        
        def retrieveAlignment(self, stateL):
            seq = self.seq 
            alignment_seq = []
            for i in range(len(seq)):
                state = stateL[i]
                if state == 'M':
                    alignment_seq.append(seq[i])
                else:
                    alignment_seq.append('-')








