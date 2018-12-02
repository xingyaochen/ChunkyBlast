import sys, fasta
import hmmmodel

import hmmAlign
    
#### Main

if __name__ == "__main__":


    if len(sys.argv) !=3:
        sys.stderr.write("""
        Usage: python seedAlign.fa protSeq.fa

""")
        sys.exit(-1)

    seedAlignFN = sys.argv[1]
    protSeqFN = sys.argv[2]

    # create profile HMM based on seed alignment
    model=hmmmodel.HmmModel(seedAlignFN)

    # load db
    sseq = fasta.load(protSeqFN)[0][1]
    alignment= hmmAlign.HmmAlign(sseq,model)
    score, L =alignment.viterbi()
    print(score)

    print(alignment.backTrack(L[0], L[1], L[2]))
    # print

    

