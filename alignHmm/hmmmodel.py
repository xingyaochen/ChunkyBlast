from fasta import load
from collections import Counter

aaList=['G','A','L','M','F','W','K','Q','E','S','P','V','I','C','Y','H','R','N','D','T']
class HmmModel:
    def __init__(self,filename):
        self.filename=filename
        self.emmit={}
        self.transit={}
        self.num_states=0
        self.info=None
        self.matchStates=None
        self.findStates()

    def findStates(self):
        """
        Load the fasta alignment file and return a list of column
        indices of the Match states
        """
        info=load(self.filename)
        self.info=info
        numCol=len(info[0][1])
        matches=[]
        
        
        for i in range(numCol):
            count=0
            for x in range(len(info)):
                if info[x][1][i]!='.':
                    count+=1
            if count>=0.5*len(info):
                matches.append(i)
        self.matchStates=matches
        return 

    def calcAllEmitProb(self):
        for matchStateNum in range(1,len(self.matchStates)+1):
            self.calcEmitProb(matchStateNum)
            
    def calcEmitProb(self,stateNum):
        """
        calculate the emmision probability for matching state with number stateNum
        fill the entry in the dictionary
        """
        # assume stateNum indexes sequnece position
        # get the corresponding column
        columnData = [seq[1][self.matchStates[stateNum-1]] for seq in self.info]

        # make a dict of aa counts 
        countsD = Counter(columnData)

        # get the number of gaps, which is excluded in the denominator
        gaps = countsD.get('.')

        #calulate denominator
        denominator = len(columnData) + 20
        if gaps:
            denominator -= gaps
        stateNumD = {}

        # put probs in a dictionary
        for aa in aaList:
            if aa in countsD:
                stateNumD[aa]=(countsD[aa]+1)/float(denominator)
            else:
                stateNumD[aa] = 1.0/denominator

        # put dictionary into meta dictionarys
        self.emmit[stateNum] = stateNumD

    def calAllTransProb(self):
        #lets encode begin as M0, end as M_len(self.matchStates)+1
        for i in range(len(self.matchStates)+1):
            self.calcTransProb(i,i+1)
      



    def calcTransProb(self, state1Num, state2Num):
        """
        calculate the transition probability for going from state1 to state2
        """
        # state: (type, pos), ex: ('M', 0), ('D', 1)
        transitionD = self.tallyTransD(state1Num, state2Num) 
        if state2Num<=len(self.matchStates):
            PossibleTransTo = [ "M" , "I", "D"]
        else:
            PossibleTransTo = [ "M" , "I"]
        
        if state1Num==0:
            PossibleTransFrom=[ "M" , "I"]
        else:
            PossibleTransFrom = [ "M" , "I", "D"]

        # calculate all denominator types
        #plus 3 for number of states to transition to
        demominatorM = sum([transitionD["M"+trans] for trans in PossibleTransFrom]) + len(PossibleTransTo)  # sum of MM, MI, MD
        if state1Num!=0:
            demominatorD = sum([transitionD["D"+trans] for trans in PossibleTransFrom]) + len(PossibleTransTo)  # sum of DM, DI, DD
        demominatorI = sum([transitionD["I"+trans] for trans in PossibleTransFrom]) + len(PossibleTransTo)  # sum of IM, II, ID
        
        # for-loop to make all possible transition types
        for i in range(len(PossibleTransFrom)):
            fromState = PossibleTransFrom[i]
            for j in range(len(PossibleTransTo)):
                toState = PossibleTransTo[j] 

                # depending on fromState, select the demominator to use
                if fromState == "M":
                    demominator = demominatorM 
                elif fromState == "I":
                    demominator = demominatorI 
                elif fromState == "D":
                    demominator = demominatorD  

                # get numerator
                numerator = transitionD[fromState+toState]

                # make keys for transit dict
                fromStateTuple = (fromState, state1Num)

                # special Insert case????
                if toState == "I":
                    toStateTuple = (toState, state1Num)
                else:
                    toStateTuple = (toState, state2Num)

                # put probs into the transition prob dictionary 
                if fromStateTuple not in self.transit:
                    self.transit[fromStateTuple] = {}
                #plus 1 in numerator for pseudocount
                self.transit[fromStateTuple][toStateTuple] = float(numerator+1)/demominator


            

    def tallyTransD(self, state1Num, state2Num):
        """
        tally up all transition types going from match state 1 to match state 2
        """
        transitionD = {}
        #special edge case for from beginning states to match1 etc.
        if state1Num==0:
            #from begin
            transitionD = {"II": 0, "IM": 0, "ID":0, "MM":0, "MI":0, "MD":0}
            nextPos = self.matchStates[state2Num-1] 
            state2Col = [seq[1][nextPos] for seq in self.info]
            if nextPos==0:
                transitionD["MD"]=state2Col.count('.')
                transitionD["MM"]=len(self.info)-state2Col.count('.')
            if nextPos!=0:
                
                initialStateInfo=[seq[1][:nextPos] for seq in self.info]

                ###TODO
                for i in range(len(self.info)):
                    #loop through each sequence
                    stateTransType = ""
                    insertion=initialStateInfo[i]
                    numInserts = len(insertion) - insertion.count('.')
                    aa2 = state2Col[i]
                    if insertion!=0:
                        stateTransType ="MI"
                    elif aa2 == ".":
                        stateTransType+="MD" 
                    elif aa2 in aaList:
                        stateTransType+="MM" 

                    transitionD[stateTransType] += 1
                    
                    if aa2 in aaList:
                        stateTransType ="IM"
                    elif aa2 == ".":
                        stateTransType = "ID" 
                    transitionD[stateTransType] += 1
                
                    stateTransType = "II"
                    if numInserts!=0:
                        transitionD[stateTransType] += numInserts-1
                
            return transitionD
        
        if state2Num==len(self.matchStates)+1:
            #special case for from last match to ending
            transitionD = {"II": 0, "IM": 0, "MM":0, "MI":0,"DI":0, "DM":0}
            curPos = self.matchStates[state1Num-1] 
            state1Col = [seq[1][curPos] for seq in self.info]
            if curPos==len(self.info[0][1])-1:
                #if the last matching column is the last column
                transitionD["DM"]=state1Col.count('.')
                
                transitionD["MM"]=len(state1Col)-state1Col.count('.')
            else:
                endingStateInfo=[seq[1][curPos+1:] for seq in self.info]
                for i in range(len(self.info)):
                    #loop through each sequence
                    stateTransType = ""
                    insertion=endingStateInfo[i]
                    numInserts = len(insertion) - insertion.count('.')
                    aa1 = state1Col[i]
                    if insertion!=0:
                        stateTransType ="MI"
                    else:
                        stateTransType+="MM" 

                    transitionD[stateTransType] += 1
                    
                    if numInserts!=0:
                        stateTransType ="IM"
                        transitionD[stateTransType] += 1
                
                    stateTransType = "II"
                    if numInserts!=0:
                        transitionD[stateTransType] += numInserts-1
            return transitionD



            

                

        curPos = self.matchStates[state1Num-1] 
        nextPos = self.matchStates[state2Num-1] 
        state1Col = [seq[1][curPos] for seq in self.info]
        state2Col = [seq[1][nextPos] for seq in self.info]
        stateMidCol = [seq[1][curPos+1:nextPos] for seq in self.info]
        transitionD = {"II": 0, "IM": 0, "ID":0, "MM":0, "MI":0, "MD":0, "DD":0, "DI":0, "DM":0}
        for i in range(len(state1Col)):
            stateTransType = ""
            aa1 = state1Col[i]
            aa2 = state2Col[i]
            if aa1 in aaList:
                stateTransType +="M"
            elif aa1 == ".":
                stateTransType+="D"
            if aa2 in aaList:
                stateTransType +="M"
            elif aa2  == ".":
                stateTransType+="D"
            if stateTransType in transitionD:
                transitionD[stateTransType] +=1
            else:
                transitionD[stateTransType] = 1 
        if stateMidCol[0]:
            for i in range(len(state1Col)):
                stateTransType = ""
                aa1 = state1Col[i]
                aa2 = state2Col[i]
                if aa1 in aaList:
                    stateTransType ="MI"
                elif aa1 == ".":
                    stateTransType+="DI" 

                transitionD[stateTransType] += 1
              
                
                if aa2 in aaList:
                    stateTransType ="IM"
                elif aa2 == ".":
                    stateTransType = "ID" 
                transitionD[stateTransType] += 1
               
                numInserts = len(stateMidCol[i]) - stateMidCol[i].count('.')
                stateTransType = "II"
                if numInserts:
                    transitionD[stateTransType] += numInserts-1
            return transitionD



    


            

