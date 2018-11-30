import sys, fasta
import hmmmodel

def alignmentScores(model,sseq):
    """
    implements the viterbi global alignment and calculate the alignment
    score for sseq against the model.
    Align sequence to model using affine gap type global alignment algorithm.
    """
    totLength=len(model.info[0][1])
    
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

    scores=alignmentScores(model,sseq)

