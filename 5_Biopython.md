# Bioinformatics x Python

## Welcome

Congratulations! After your hard work honing your Unix and Python skills, you've been selected for a bioinformatic research project! Soon, you're going to be trusted with genetic data from an unknown source. Your goals will be to use bioinformatics to 1) analyze features of the data and 2) determine where this data came from.

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

Of course not. Let's find out [how things are really done](6_BiopythonV2.md).
