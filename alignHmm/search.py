import sys, fasta
import hmmmodel

import hmmAlign


#### Main

if __name__ == "__main__":


    if len(sys.argv) !=4:
        sys.stderr.write("""
        Usage: python seedAlign.fa protDB.fa outfile.txt

""")
        sys.exit(-1)

    seedAlignFN = sys.argv[1]
    dbSeqFN = sys.argv[2]
    outFN = sys.argv[3]

    # create profile HMM based on seed alignment
    seedAlignL = [seq for header,seq in fasta.load(seedAlignFN)]

    # load db
    dbL = fasta.load(dbSeqFN)

    outL=[]
    model=hmmmodel.HmmModel(seedAlignFN)
    for hd,sseq in dbL:
        # load db
        print(hd)
        alignment= hmmAlign.HmmAlign(sseq,model)

        score=alignment.subtractShuffleMean()
        
        outL.append((score,hd))

    outL.sort(reverse=True) # sort by score, high to low

    # write sorted output to file
    f=open(outFN,'w')
    for sc,hd in outL:
        print(sc,hd,file=f)
    f.close()
    
