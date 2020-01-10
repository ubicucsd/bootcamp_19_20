# Alignment

## EC2: ec2-18-191-106-205.us-east-2.compute.amazonaws.com

#### Skill: More python and some knowledge of alignment methodology

Every bioinformatics tool needs to start somewhere, and beyond quality control the first step most tools require is to make sequences comparable by aligning them. Alignment basically just put regions that are similar to each other closer to each other, and inserts gaps (denoted by '-') where there may have been an insertion/deletion event. Identifying regions of similarity can help us identify structural, functional, or evolutionary relationships between sequences.  

There are probably hundreds of alignment methods specialized for different types of information (DNA vs amino acid), different uses (finding common domans vs comparing homologous genes), and different special cases (antibody sequences, viral sequences). If you are interested in learning details about what a few of the methods that exist and how they work (this is BENG 181 material): the most generic DNA alignment method is a pairwise dynamic programming method called [Smith Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm), alignment of many sequence (multiple sequence alignment) is done by [Fast Fourier Transforms](https://en.wikipedia.org/wiki/MAFFT), and alignment for accurate amino acid sequence homology is done by [hmm profile alignment](http://www.biology.wustl.edu/gcg/hmmanalysis.html), 

**This lesson will focus on hands on learning of three types of alignment: global, local, and multiple sequence alignments.**

## Why does an alignment look like it does?

In order to see why alignments are useful on a larger scale, let us start by looking at them from a smaller scale. 

### Global Alignment

In this case, the word "global" just means that the entire first string is aligned as best as possible to the entire second string. 

First, make a file called ```aligners.py``` and import the ```pairwise2``` module from the BioPython library and the ```format_alignment``` method from the ```pairwise2``` module. **HINT: look at part 4 to find the correct import syntax**

Now, define two small strings containing DNA sequences. I recommend using ```TGCCTTAG``` and ```TGCTTGC``` for an easy to look at example. Call the default pairwise alignment method, called 

```pairwise2.align.globalxx()``` 

on those two strings. The alignment method returns a list of the most high scoring(good) alignments. You can print those out by iterating through the alignments with a for-each loop and print them out one by one. In order to make the alignments look a little nicer, you can put them through a formatting method before printing:

```print(format_alignment(*alignment_name))```

What is that asteriks? In different languages an asteriks has different meanings, but in Python it denotes a variable quantity of arguments. In other words, methods that want this asteriks are capable of accepting either one, maybe two, or maybe more arguments(inputs). 

### What exactly is a good score? Why is one alignment better than another? 

Glad you asked! If you look at the score printed out below each alignment, you will notice that the score is coincidentially identical to the number of nucleotides that match between the two aligned sequences. The way we determine the alignment with the optimal score is beyond the scope of what I plan to discuss, but it is good to know that alignments are mostly determined by the way a program decides on alignment scores. The alignment we tried gives gaps, insertions, and deletions a score of zero and matches a score of one. Other alignments may have more complex behaviors, including negative scores for insertions/deletions, custom scores based on which nucleotide is mismatched with which, and more. There are several modes of alignment available on BioPython, each with customizable scoring. Look [over here](http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html) if you want to see the range of what BioPython has to offer.

### Local Alignment

The word "local" means the best alignment between two *subsequences*. This method can be used when looking for a conserved gene between two otherwise very different organisms. Aligning the messy differences between two different species is not useful, but finding two subsequences (two genes) that the species have in common without aligning the whole sequences can be very useful. 

Open up ```aligners.py``` and create strings ```ATGCGGCCATTATAAGCGGTCTCT``` and ```GATTATAATT```. Now, let's compare how these strings align with globally vs locally. The point is best illustrated when gaps and mismatches are penalized, so indicate the scoring system of the alignment like this: 

```pairwise2.align.localms(str1, str2, 1, -1, -1, -1)```
```pairwise2.align.globalms(str1, str2, 1, -1, -1, -1)```

The numbers at the end indicate +1 for matches, -1 for mismatches, -1 for opening a gap, and -1 for extending a gap. 

Looking at the matching nucleotides in global and local alignments of these string, which one makes more sense? Does it make sense to care about the matches the global alignment has at the very end of the sequences? 

## Multiple Sequence alignment 

Multiple sequence alignment aligns multiple sequences, but its inner workings are bit complicated (my way of saying I do not know them well enough to teach them) so we are just going to look at them from a distance. This type of alignment is used on a large number of more or less related sequences in order to infer homology and build evolutionary trees. My multiple sequence aligner of choice is mafft, which we will be using in the challenge below. 
