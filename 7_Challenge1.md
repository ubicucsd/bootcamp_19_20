# Challenge 1

### EC2: ec2-13-58-73-19.us-east-2.compute.amazonaws.com

## Welcome
Welcome to Week 4 of the Crash Course! This will be our final session this quarter; good luck with finals, and we'll see you in Winter.

This week's lesson is, in some ways, a one-question "midterm" that covers the skills you've learned in Weeks 1-3. It is designed to be more representative of the difficulty level of problems you'd be asked to solve in real world applications rather than a tutorial explaining how things are done. It is more challenging than anything we've done so far, so don't be discouraged if you have to think things over (and fail) before getting things right!

## Instructions
### Coverage
We assume that you've already covered (and only covered) all three of our previous sessions. For your reference, they're linked here: [Week 1](/1_Welcome.md) [Week 2](/3_AdvancedTerminal.md) [Week 3](/5_Biopython.md). The techniques for solving the problem below are, with 100% certainty, covered in the previous three lessons. 

*Note: This does not mean that they are explicitly mentioned in one of our previous lessons. Recall that we emphasize finding specific commands, options, and tools on your own. However, we have covered the core concepts before.*

### Logistics
The problem below does not specify *how* you should solve it. Instead, we'll present you with a task to perform and you'll need to determine the output.

Additionally, you'll never need to download any software. Everything you need is on EC2.

### Getting Help
In the real world, you're often on your own when you're stuck and need help. Reading through Google searches, documentations, tutorials, and manuals is a real life skill. 

However, we don't just want to throw you off the deep end. If you do get stuck on something and you can't figure it out, you can draft and email and send it to Sabeel. Provide a specific description of 1) what you're trying to accomplish 2) what issue you're running into 3) what you've already tried to troubleshoot it. If we find many common issues, we'll send out an FAQ/hints email on Sunday.

### Good luck, have fun!

## The Problem
We've previously described how and why we sequence DNA. The raw output of these sequencing machines are [FASTQ](https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html) files. Recall that each base that a FASTQ file contains is paired with a quality score that's encoded in [ASCII](http://www.asciitable.com/) (proportional to how confident the machine is that the base was identified correctly). We, of course, do not want to have our downstream analyses tainted by low-confidence sequences. 

There are various techniques for cleaning up such raw data. Your goal here is to implement one possible method: We will only keep reads **if the average base quality is (strictly) greater than 68**.

As an example, imagine our FASTQ file contains only two sequences, each with three characters:
```
@SAMPLE1
AGG
+
DEF
@SAMPLE2
TGC
+
BCD
```

@SAMPLE1 has an average quality score of 69 (from the average ASCII of `DEF`). We will keep this sequence.
@SAMPLE2 has an average quality score of 67 (from the average ASCII of `BCD`). We will not keep this sequence.

The output file, therefore, would contain four lines:
```
@SAMPLE1
AGG
+
DEF
```

Make a copy of your input FASTQ file, available at `~/../smansuri/raw.fastq`. Have fun!

## Checking Your Solution
You'll be self-checking your work this week. If your output file has 3828 lines, you're done! Let us know you're done by sending Sabeel an email titled **[UBIC] Completed Week 4** containing the last line of your output file (i.e. line 3283).

## Congratulations! You've completed Week 4 AND Quarter 1 of the Bioinformatics Crash Course.
