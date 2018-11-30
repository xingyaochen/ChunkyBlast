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
        #matrices have index (j,i) where j range from 0 to len(states)+1, i ranges from 0 to len(seq)-1
        v_m, v_i, v_d = [[None]*(len(seq)+1) for i in range(len(states)+2)], [[None]*len(seq) for i in range(len(states)+2)], [[None]*len(seq) for i in range(len(states)+2)]
        #initialize:
        #matching state matrix
        v_m[0][0]={"prob": 0, "prev": None}
        # for stNum in range(len(states)+1):
        #     v_m[stNum][0] = {"prob": 0, "prev": None}
        #     v_i[stNum][0] = {"prob": 0, "prev": None}
        #     v_d[stNum][0] = {"prob": 0, "prev": None}
        qxi = 666 
        for i in range(1, len(seq)-1):
            aa = seq[i]
            for j in range(1, len(states)):
                e_m = hmmmodel.getEmitProb(st, aa) 

                amm = hmmmodel.getTransitProb(('M', i), ('M', i+1))
                aim = hmmmodel.getTransitProb(('I', i), ('M', i))
                adm = hmmmodel.getTransitProb(('D', i), ('M', i+1))

                ami = hmmmodel.getTransitProb(('M', i), ('I', i))
                aii = hmmmodel.getTransitProb(('I', i), ('I', i))
                adi = hmmmodel.getTransitProb(('D', i), ('I', i))

                amd = hmmmodel.getTransitProb(('M', i), ('D', i+1))
                aid = hmmmodel.getTransitProb(('I', i), ('D', i))
                add = hmmmodel.getTransitProb(('D', i), ('D', i+1))

                # for just v_m 

                m_threeProbs = [v_m[j-1][i-1]["prob"] + amm , v_i[j-1][i-1]["prob"] + aim, v_d[j-1][i-1]["prob"] + adm] 

                max_prob_m = np.amax(m_threeProbs)
                prev_state_m = (self.states[np.argmax(m_threeProbs)], i-1, j-1)
                
                vm[i][j] = {"prob": math.log(e_m / qxi)  + max_prob_m, "prev": prev_state_m}


                # TODO: do the same for v_i and v_d

            # TODO: take care of edge cases
            
        bestScore=v_m[len(states)+1][len(seq)-1]
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








