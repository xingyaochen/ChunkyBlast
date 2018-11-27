from fasta import load
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
                matches+=[i]
        self.matchStates=matches
        return 

    def calcEmitProb(self,stateNum):
        """
        calculate the emmision probability for matching state with number stateNum
        fill the entry in the dictionary
        """
    
    def calcTransProb(self,state1,state2 ):
        """
        calculate the transition probability for going from state1 to state2
        """
    


            

