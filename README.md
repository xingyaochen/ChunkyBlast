# ChunkyBlast

## Materials
The file tim-PF00121_seed.fa is a seed alignment for the triosephosphate isomerase (TIM) protein family from PFAM. You will build a profile hidden Markov model from this alignment.

(A note for the curious: triosephosphate isomerase is a glycolytic enzyme, widely distributed across the tree of life.  This particular enzyme is very effective, accelerating its reaction 10^9 fold. It has been termed a "perfect catalyst" in the sense that the rate of this reaction appears to be limited by the speed at which the product diffuses away from the enzyme, rather than the speed of activity of the enzyme itself.)

The file tim-searchSeq.fa represents a triosephosphate isomerase sequence not in the seed alignment. It can be used for testing the alignment algorithms you will write.

The multi fasta file unknownProts.fa is a set of unknown protein sequences which you can search with your model, trying to find those with significant similarity to the TIM family

The python files alignOne.py and search.py will provide the interface for aligning one sequence against the model, or searching a whole set. You will want to modify them, and also create several other support files with functions and classes which they import. Once your code is complete, they can be used like this:

python alignOne.py tim-PF00121_seed.fa tim-searchSeq.fa
python search.py tim-PF00121_seed.fa unknownProts.fa timOut.txt

The file unknownProts.fa contains 101 sequences. A couple of these are actually related to triosephosphate isomerase. Because we're not doing any sort of statistics with our searches, you won't have any way of knowing what log odds score is high enough to definitively indicate homology. However if things are working right, there should be several sequences at the top of your list with much higher scores than the rest.

## The Model

Your first task is to create a probabilistic model of a protein family based on a seed  alignment. In class we discussed how to get probabilities out of a seed alignment. Chapter 5 of Durbin et al. (linked above) also contains a discussion.

To implement this, you will need to think about how best to represent the many different probabilities included in a profile HMM. You will then need to create functions to read through a seed alignment calculating the probabilities and filling up your data structure. You will also need to create functions to get various probabilities out of the data structure when they are needed (ie by an alignment algorithm).

One way to do this is to create a python class to represent your profile HMM. Internally the class will use things like dictionaries and lists to hold the probabilities. But it will also have methods built into it that allow users to extract those probabilities simply (and without needing to know the details of how they are represented).

## Aligning a sequence to the model

The search application requires us to obtain a score for aligning a sequence against a model.

Viterbi and backtracking. The most basic way to get a score is by implementing the Viterbi algorithm (global alignment) for profile HMMs. Strictly speaking, the search application doesn't require us to do backtracking. However being able to look at the alignments directly is extremely useful for debugging purposes, and you will almost certainly want to do backtracking.

Getting a log odds score. You will definitely want to be working in log space. So the scores that you are dealing with in the alignment algorithm will be log probabilities, and the final score from Viterbi will be log of the probability of the best path. This kind of a score is sensitive to sequence length (the probability goes down the longer the sequence), which is bad for the search application. A simple solution is to calculate a second score for a random model. The easiest way to do this is to (repeatedly) shuffle the sequence you are aligning, getting the average score of the shuffled sequence against the model. Subtracting this from the real (unshuffled) score gives a log odds, which is not sensitive to length and is suitable for the search application.

On the forward algorithm. For this project, it is sufficient to implement Viterbi and use that in search.py. However, if you are so motivated, you could additionally implement the forward algorithm which considers all possible alignment paths, and is very powerful for finding distant protein similarities. As was the case for Viterbi, you will want to work in log space here. This makes the forward algorithm a bit tricky. In it you are required to add probabilities. Doing that directly would mean coming out of log space, which is not a good idea (either for accuracy or efficiency). Fortunately, numpy has a function, numpy.logaddexp which takes two numbers in log, and does the equivalent of 1) exponentiating them, 2) adding the result, and 3) putting it back in log.

## Submitting your assignment

Once you are done, please write me a text file, notes.txt, giving an overview of how you set up your code, and what files it's located in. Please also discuss the strategy you took for testing your code, and any other comments that are relevant. Put notes.txt in your working directory, and then compress the whole thing back into a zipfile (also call this hmmAlign.zip). Submit the zip file at the submission site.


## Further reading
Albery, W.J. and Knowles, J.R., 1976. Evolution of enzyme function and the development of catalytic efficiency. Biochemistry, 15(25), pp.5631-5640.

Durbin, R., Eddy, S.R., Krogh, A. and Mitchison, G., 1998. Biological sequence analysis: probabilistic models of proteins and nucleic acids. Cambridge university press.

Nagano, N., Orengo, C.A. and Thornton, J.M., 2002. One fold with many functions: the evolutionary relationships between TIM barrel families based on their sequences, structures and functions. Journal of molecular biology, 321(5), pp.741-765.
