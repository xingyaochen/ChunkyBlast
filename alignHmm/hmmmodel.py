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


    def calcTransProb(self, state1Num, state2Num):
        """
        calculate the transition probability for going from state1 to state2
        """
        # state: (type, pos), ex: ('M', 0), ('D', 1)
        transitionD = self.tallyTransD(state1Num, state2Num) 

        allPossibleTrans = [ "M" , "I", "D"]
        # calculate all denominator types
        demominatorM = sum([transitionD["M"+trans] for trans in allPossibleTrans]) + 3  # sum of MM, MI, MD
        demominatorD = sum([transitionD["D"+trans] for trans in allPossibleTrans]) + 3  # sum of DM, DI, DD
        demominatorI = sum([transitionD["I"+trans] for trans in allPossibleTrans]) + 3  # sum of IM, II, ID
        
        # for-loop to make all possible transition types
        for i in range(len(allPossibleTrans)):
            fromState = allPossibleTrans[i]
            for j in range(len(allPossibleTrans)):
                toState = allPossibleTrans[j] 

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
                self.transit[fromStateTuple][toStateTuple] = float(numerator)/demominator


            

    def tallyTransD(self, state1Num, state2Num):
        """
        tally up all transition types going from match state 1 to match state 2
        """
        transitionD = {}
        curPos = self.matchStates[state1Num-1] 
        nextPos = self.matchStates[state2Num-1] 
        state1Col = [seq[1][curPos] for seq in self.info]
        state2Col = [seq[1][nextPos] for seq in self.info]
        stateMidCol = [seq[1][curPos+1:nextPos] for seq in self.info]
        transitionD = {"II": 1, "IM": 1, "ID":1, "MM":1, "MI":1, "MD":1, "DD":1, "DI":1, "DM":1}
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


# we may not need all this ...
    def transProbMtoM(self, curPos, nextPos):
        pass 
    
    def transProbMtoI(self, curPos, nextPos):
        pass 
    
    def transProbMtoD(self, curPos, nextPos):
        pass 
    
    def transProbDtoM(self, curPos, nextPos):
        pass 
    
    def transProbDtoI(self, curPos, nextPos):
        pass 
    
    def transProbDtoD(self, curPos, nextPos):
        pass 
    
    def transProbItoM(self, curPos, nextPos):
        pass 
    
    def transProbItoI(self, curPos, nextPos):
        pass 
    
    def transProbItoD(self, curPos, nextPos):
        pass 

    


            

