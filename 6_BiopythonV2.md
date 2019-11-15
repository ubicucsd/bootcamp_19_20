# Bioinformatics x Python (v2)

## Welcome

Congratulations! After your hard work honing your Unix and Python skills, you've been selected for a bioinformatic research project! Soon, you're going to be trusted with genetic data from an unknown source. Your goal will be to use bioinformatics to determine where this data came from.

## Warm Up (Packages)

Before you're allowed to work with real genetic data, you'll need to familiarize yourself with Python packages. Think of packages as Python programs someone else has written that you can borrow and use for your own purposes.

For example, earlier you determined the GC content of your sequence by writing your own algorithm. Instead, it's a lot easier to find a Python package that already has a built-in GC counter, and just run it by doing something like `sequence.count_gc()`, for example.

The Python package we'll be using for our bioinformatics purposes is called Biopython. This is already installed on our cluster. If it wasn't, you could install it by doing (DO NOT ACTUALLY DO THIS): `pip install --user biopython`.

Instead, prove to yourself that Biopython is already on the cluster by doing:
```
pip list
```

## Warm Up (Biopython)

### Seq Objects

Much like the fundamental building block of genetics is DNA, the fundamental building block of Biopython (at least for what we're doing) is the Seq object. 

Great, we have Biopython. Let's figure out what it does and how to use it.

## The First Glimpse

### Access

You'll be using the same genetic data as before. This time, instead of having various files for various methods, we're going to be using just one file. Isn't that handy? Copy it to your directory:
```
cp ~/../smansuri/biopython.py .
```

### File Extenstion

Confirm that you renamed the initisl file as *unidentified.fasta*. If you did not, review [FASTA format](https://www.genomatix.de/online_help/help/sequence_formats.html#FASTA)

## Transcription Simulation

### Reverse Complement

Before, you had to write the logic for the reverse complement yourself. Now, just use the Seq object's `reverse_complement()` method.

### RNA Transcription

Similar to above, use the `transcribe()` method on the reverse complement. Voila, one line and you're done.

## A Potential Breakthrough

Recall the two markers that could suggest your unknown sample may, in fact, be related to the one seen in rats:
1. The sequence "ATGGAGCTGACTGTGGAGGCATG" is often present.
2. The GC Content is above 55%.

Use the Seq object's `find()` method for #1, and explore the Biopython [documentation](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf) to find a method for #2.

## BLAST



## Congratulations! You've completed Week 3 of the Bioinformatics Crash Course.

## Credits
Exercises are adapted from [Rosalind](http://rosalind.info), the official [Biopython tutorial textbook](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf), and a [course](http://disi.unitn.it/~teso/courses/sciprog/python_biopython_exercises.html) from the University of Trento.












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


