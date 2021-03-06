import math
import fasta
import numpy as np
import numpy as np
from random import shuffle



class HmmAlign:
    def __init__(self, sseq, hmmmodel):
        self.seq = sseq
        self.hmmmodel = hmmmodel 
        self.states = ['M', 'I', 'D']

    def viterbi(self):
        #searching seq
        seq = self.seq
        states = self.hmmmodel.matchStates
        #matrices have j rows ,i columns where j range from 0 to len(states)+1, i ranges from 0 to len(seq)-1
        v_m, v_i, v_d = [[{}]*(len(seq)+1) for i in range(len(states)+2)], [[{}]*(len(seq)+1) for i in range(len(states)+2)], [[{}]*(len(seq)+1) for i in range(len(states)+2)]
        #initialize:
        #initialize M matrix
        # transition from I0 is allowed
        v_m[0][0]={"prob": 0, "prev": None}
        for i in range(1,len(seq)+1):
            #begin state M0 cannot match to sequence nucleotide positions 
            v_m[0][i]= {"prob": float('-inf'), "prev": None}
        for j in range(1,len(states)+2):
            #can't match fake position to any matching state
            v_m[j][0]= {"prob": float('-inf'), "prev": None}
        #end initializing M matrix

        #initialize I matrix: 
        # transition from I0 is allowed
        v_i[0][0]={"prob": 0, "prev": None}
        for j in range(1, len(states)+2):
            #can't match fake position to any insertion state
            v_i[j][0]=  {"prob": float('-inf'), "prev": None}
        #end initialize I matrix
        #initialize D matrix:
        for i in range(1,len(seq)+1):
            #D0 doesn't exit
            v_d[0][i]= {"prob": float('-inf'), "prev": None}
        for j in range(len(states)+2):
            #can't match fake position to any insertion state
            v_d[j][0]=  {"prob": float('-inf'), "prev": None}
        #end initializing D matrix
        for i in range(1, len(seq)+1): # i indexes the amino acid
            aa = seq[i-1]
            for j in range(len(states)+2): # j indexes the state
                if j == 0: # fill up only v_i 
                    e_i = self.hmmmodel.getEmitProb(("I", j), aa) 
                    ami = self.hmmmodel.getTransitProb(('M', j), ('I', j))
                    aii = self.hmmmodel.getTransitProb(('I', j), ('I', j))
                    adi = self.hmmmodel.getTransitProb(('D', j), ('I', j))
                    i_threeProbs =  [v_m[j][i-1]['prob'] + ami , v_i[j][i-1]['prob'] + aii, v_d[j][i-1]['prob'] + adi]
                    max_prob_i = np.amax(i_threeProbs) + e_i
                    prev_state_i = (self.states[np.argmax(i_threeProbs)], j, i-1)
                    v_i[j][i] = {"prob": max_prob_i, "prev": prev_state_i}
                else: # fill up all dp tables
                    if j == len(states)+1:
                        e_m = 0
                    else:
                        # get all log probabilities
                        # emissions
                        e_m = self.hmmmodel.getEmitProb(('M', j), aa)
                        e_i = self.hmmmodel.getEmitProb(("I", j), aa)  
                        # transitions
                        ami = self.hmmmodel.getTransitProb(('M', j), ('I', j))
                        aii = self.hmmmodel.getTransitProb(('I', j), ('I', j))
                        adi = self.hmmmodel.getTransitProb(('D', j), ('I', j))
                        amd = self.hmmmodel.getTransitProb(('M', j-1), ('D', j))
                        aid = self.hmmmodel.getTransitProb(('I', j-1), ('D', j))
                        add = self.hmmmodel.getTransitProb(('D', j-1), ('D', j))

                        # find the max for each dp table
                        i_threeProbs = [v_m[j][i-1]['prob'] + ami , v_i[j][i-1]['prob'] + aii, v_d[j][i-1]['prob'] + adi]
                        d_threeProbs = [v_m[j-1][i]['prob'] + amd , v_i[j-1][i]['prob'] + aid, v_d[j-1][i]['prob'] + add]
                        max_prob_i = np.amax(i_threeProbs) + e_i
                        max_prob_d = np.amax(d_threeProbs)

                        # get info for backtracking
                        prev_state_i = (self.states[np.argmax(i_threeProbs)], j, i-1)
                        prev_state_d = (self.states[np.argmax(d_threeProbs)], j-1, i)
                        
                        # input into dp table
                        v_i[j][i] = {"prob": max_prob_i, "prev": prev_state_i}
                        v_d[j][i] = {"prob": max_prob_d, "prev": prev_state_d}


                    # for v_m dp table, do the same
                    amm = self.hmmmodel.getTransitProb(('M', j-1), ('M', j))
                    aim = self.hmmmodel.getTransitProb(('I', j-1), ('M', j))
                    adm = self.hmmmodel.getTransitProb(('D', j-1), ('M', j))
                    m_threeProbs = [v_m[j-1][i-1]['prob'] + amm , v_i[j-1][i-1]['prob'] + aim, v_d[j-1][i-1]['prob'] + adm] 
                    max_prob_m = np.amax(m_threeProbs)
                    prev_state_m = (self.states[np.argmax(m_threeProbs)], j-1, i-1)
                    v_m[j][i] = {"prob": e_m + max_prob_m, "prev": prev_state_m}

        # best alignment score is the last match state
        bestScore=v_m[len(states)+1][len(seq)]['prob']  
        return bestScore, [v_m, v_i, v_d]

    def subtractShuffleMean(self):
        """
        viterbi alignment score substract by the randomized score
        for normalization
        """
        score,_=self.viterbi()
        score_L=[]
        l=list(self.seq)
        for i in range(20):
            # shuffle the sequence to randomize
            shuffle(l)
            alignment1= HmmAlign(l,self.hmmmodel)
            random_score, _ = alignment1.viterbi()
            score_L.append(random_score)
        meanScore=np.mean(score_L)
        return score-meanScore

    def backTrack(self, v_m, v_i, v_d):
        """
        Backtack to visualize the alignment result
        """
        seq = self.seq
        states = self.hmmmodel.matchStates
        state, j, i = v_m[len(states)+1][len(seq)]["prev"]
        stateL = [state]
        while j >= 1 and i >= 1:
            stateL.append(state)
            if state == "M":
                state, j, i = v_m[j][i]["prev"]
            elif state == "I":
                state, j, i = v_i[j][i]["prev"]
            elif state == "D":
                state, j, i = v_d[j][i]["prev"]
            
        stateL = stateL[::-1]
        return "".join(self.retrieveAlignment(stateL))

    def retrieveAlignment(self, stateL):
        """
        Helper function to retrieve alignment given a list 
        of state transitions
        """
        seq = list(self.seq)
        alignment_seq = []
        for i in range(len(stateL)):
            state = stateL[i]
            if state == 'M' or state == 'I' :
                alignment_seq.append(seq.pop(0))
            else:
                alignment_seq.append('-')
        return alignment_seq

