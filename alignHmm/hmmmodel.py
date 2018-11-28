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


    def calcTransProb(self, state1, state2):
        """
        calculate the transition probability for going from state1 to state2
        """
        # state: (type, pos), ex: ('M', 0), ('D', 1)
        curPos = state1[1]
        nextPos = state2[1]
        state1Col = [seq[1][curPos] for seq in self.info]
        state2Col = [seq[1][nextPos] for seq in self.info]

        transitionD = {}
        # for D and M cases
        if curPos!=nextPos:
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
        # for I cases????
        return


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

    


            

