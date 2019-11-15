# Bioinformatics x Python

## Welcome

Congratulations! After your hard work honing your Unix and Python skills, you've been selected for a bioinformatic research project! Soon, you're going to be trusted with genetic data from an unknown source. Your goal will be to use bioinformatics to determine where this data came from.

## Warm Up

Before you're allowed to work with real genetic data, you'll need to demonstrate that you know [command line](/2_LinuxTerminal.md) and [Python](/4_Python.md) basics. Complete the following tasks:

1. Run the following command: `TODO`
2. Execute the Python file and follow the instructions
3. Use a command line tool called `wget` to download genetic data from [here](TODO)
4. Determine the number of lines, number of bases, and the GC content

*Confirm your progress with an instructor before moving forward.*

## The First Glimpse

### Access

Great! You've now been granted access to the data you're tasked with identifying. Download it by running the following commands:
```shell
TODO
```

As you've always done, explore the data to get a sense of what you're working with. There are no requirements of what you need to do.

### File Extenstion

The file you downloaded is called `unidentified.unknown`. Having worked with genetic data files before, you know what type of file this is. Change the file extension from `.unknown` to the appropriate extension. *For example, if you thought this was a PDF file, you would rename the file to `unidentified.pdf`.*

## RNA Simulation

It's clear that the genetic data in this file is DNA. 

Create a file that takes in the DNA data (A/T/C/G) from "test.fasta" and prints out the file as RNA data (A/U/C/G). (Note: This is not how real transcription works. Just use this simple replacement as an example).


## A Potential Breakthrough

The transcribed RNA you submitted has a sequence that's fairly similar to --TODO--. Specifically, there are two markers that would suggest your unknown sample is --:
  1. The sequence "TODO" is often present.
  2. The GC Content is between TODO.

Determine if this sequence appears in your sample, and see if the GC content is in the expected range.

## BLAST

Your colleague notices what you've been trying to do. She suggests you use an online tool called [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/BlastAlign.cgi) to compare your sequence to a database of sequences online. Use the Nucleotide Blast tool to determine what this sample is.

## Let's Try That Again

What you just did--exploring an unknown sequence--is one of the simplest and most common bioinformatic tasks. Ask yourself this: does everyone who wants to do this have to write their own algorithms, manually parse through sequences, and copy-paste into a web browser? 

Of course not. Let's find out [how things are really done](TODO).












# Bioinformatics x Python

## Welcome

Congratulations! After your hard work honing your Unix and Python skills, you've been selected for a bioinformatic research project! Soon, you're going to be trusted with genetic data from an unknown source. Your goal will be to use bioinformatics to determine where this data came from.

## Warm Up

Before you're allowed to work with real genetic data, you'll need to demonstrate that you know [command line](/2_LinuxTerminal.md) and [Python](/4_Python.md) basics. Complete the following tasks:

1. Run the following command: `TODO`
2. Execute the Python file and follow the instructions
3. Use a command line tool called `wget` to download genetic data from [here](TODO)
4. Determine the number of lines, number of bases, and the GC content

*Confirm your progress with an instructor before moving forward.*

## The First Glimpse

### Access

Great! You've now been granted access to the data you're tasked with identifying. Download it by running the following commands:
```shell
TODO
```

As you've always done, explore the data to get a sense of what you're working with. There are no requirements of what you need to do.

### File Extenstion

The file you downloaded is called `unidentified.unknown`. Having worked with genetic data files before, you know what type of file this is. Change the file extension from `.unknown` to the appropriate extension. *For example, if you thought this was a PDF file, you would rename the file to `unidentified.pdf`.*

## Transcription

It's clear that the genetic data in this file is DNA. 

Create a file that takes in the DNA data (A/T/C/G) from "test.fasta" and prints out the file as RNA data (A/U/C/G). (Note: This is not how real transcription works. Just use this simple replacement as an example).

## Let's Try That Again

What you just did--determining the source of an unknown sequence--is one of the simplest and most common bioinformatic tasks. Ask yourself this: does everyone who wants to do this have to write their own algorithms, manually parse through sequences, and copy-paste into a web browser? 

Of course not. Let's find out [how things are really done](TODO).




## Packages

An important feature of Python is its packages. Think of packages as Python programs someone else has written that you can borrow and use for your own purposes.

For example, you may want to take the first sequence of "test.fasta" and use NCBI BLAST (learn more about the BLAST database [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi)) to determine what organism the sequence may have come from. You wouldn't want to write your own code to connect to the database... instead you can use a nifty package to do that for you. Let's try:

First, extract the first two lines of "test.fasta" into a new file, "small.fasta":

```shell
head -n2 test.fasta >> small.fasta
```

Now, install the correct package. The one you want to use is called [BioPython](https://biopython.org/). Use the Python program installer (called pip) to install the package to your user:

```shell
pip install --user biopython
```

Great! BioPython is now available on your account. Let's use it. Create a new file "Blast.py", and add the following code:

```python
from Bio.Blast import NCBIWWW

fasta_string = open("small.fasta").read()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string, format_type='Text', hitlist_size=1)
print result_handle.read()
```

Let's walk through what this does. The first line takes the programs we want for NCBI from BioPython and prepares them to be used. The second line reads in the "small.fasta" file. The third line is the most important: it takes the genetic data, connects to the NCBI BLAST database, searches for matches, and then returns the result from the database. Finally, the last line prints that result.

You can try running "Blast.py" now if you'd like, but I'd recommend coming back to this after completing the next 3 exercises.




**Make sure you don't attempt to transcribe the header lines!**

#### SOLUTION
```python
file = open("test.fasta", "r")

for line in file:
  if line[0] != ">":
    ans = ""
    for char in line:
      if char == "T":
        ans += "U"
      else:
        ans += char
    print ans
  else:
    print line
```


## Introduction

This lesson is a continuation of our previous lesson, Intro to Python, . If you're having trouble following this tutorial, are new to Python, or just want a refresher, check that lesson out first.

## Setup

We're going to start by using a Python package called [Biopython](https://biopython.org/) to perform a few common bioinformatic tasks. Biopython is already installed on the cluster. You can prove this to yourself by typing the following on the command line:

```shell
pip list
```

## Sequences

### Seq objects

The first thing you think of when you hear "genetic data" is probably a sequence of DNA, which is a string of nucleotides (A's C's G's and T's). Let's use Biopython to work with sequences. Create a file to work in.

First, import the Seq package:

```python
 from Bio.Seq import Seq
```

Next, create a Seq object (think of that as a sequence) and assign it to a dna variable:
```python
dna = Seq("AGTACACTGGT")
```

Now, print out that variable:
```python
print dna
```

### Loading sequences

In reality, you're probably not going to be typing in your A-G-C-T's--you'll be using a fasta file containing them. 

Now to parse, we can load a fasta file to work with. Create a new file and import the SeqIO package:
 
 ```python
from Bio import SeqIO
for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
   print(seq_record.id)
   print(repr(seq_record.seq))
   print(len(seq_record))
```

 Before you run this program, take a guess what will print (confirm with a neighbor!). Were you right? If not, what printed and why?
 

## Complements and Transcription

You've probably heard of complementary strands and transcription. If you haven't, here are two excellent sources: [transcription](https://www.khanacademy.org/science/biology/gene-expression-central-dogma/transcription-of-dna-into-rna/a/overview-of-transcription) and [complements](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)). Biopython can perform these actions too!

### Complements

Create a new file and import the Seq package again. Create a Seq object with a DNA strand of your choosing (A's C's G's and T's only!). Call the `complement()` method on your Seq object. What do you expect will print out? Does it?

### Transcription

**Make sure you understand transcription before moving on**

Copy and paste the following into a new file:

```python
from Bio.Seq import Seq
template_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```
You've been given a template DNA strand. Can you turn this into the transcribed RNA strand? 

**Hint: You'll need the `reverse_complement()` before you `transcribe()`**

## NCBI BLAST

We connected to the BLAST database last time (take a look if you haven't already!). The output, however, wasn't very nice. Biopython can clear that up for us. 

Let's prepare to BLAST. Instead of explicitly specifying a fasta sequence, we're going to instead use a [gi number](https://www.ncbi.nlm.nih.gov/genbank/sequenceids/). Don't worry too much about this--just know we're running BLAST on some sequence, even though we haven't specified it with A's G's C's and T's.

```python
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
blast_record = NCBIXML.read(result_handle)
for alignment in blast_record.alignments:
  for hsp in alignment.hsps:
    if hsp.expect < 8.02643e-114:
      print("****Alignment****")
      print("sequence:", alignment.title)
      print("length:", alignment.length)
      print("e value:", hsp.expect)
      print(hsp.query[0:75] + "...")
      print(hsp.match[0:75] + "...")
      print(hsp.sbjct[0:75] + "...")
```
Look through the code, and look through [this documentation](http://biopython.org/DIST/docs/api/Bio.Blast.Record.HSP-class.html) (side note: a lot of coding/bioinformatics is looking through documentation--this is a valuable exercise that you shouldn't just skip!). What do you expect will print? Run the program to check.

## Your turn!

There's only one (fairly easy) challenge for this half of today's lesson, and it comes with starter code! Your task is to write a program that asks for two DNA sequences, concatenates them, and prints out the complement and reverse complement.

Starter code:
```python
seq1 = raw_input("First sequence: ")
seq2 = raw_input("Second sequence: ")

print "Complement of strings: "
print "Reverse complement of strings: "
```

## Congratulations! You've completed Week 3 of the Bioinformatics Crash Course..

## Credits
Exercises are adapted from [Rosalind](http://rosalind.info), the official [Biopython tutorial textbook](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf), and a [course](http://disi.unitn.it/~teso/courses/sciprog/python_biopython_exercises.html) from the University of Trento.
