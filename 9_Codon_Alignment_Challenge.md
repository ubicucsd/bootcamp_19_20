# Codon Alignment Challenge

One thing biologists care a lot about is the way amino acids change through time. This is found by sequencing DNA at several timepoints, codon aligning various timepoints, and comparing timepoints to see which amino acids change at what time. Why don't we just convert to amino acids and align there? Because different reading frames and deletions cause suboptimal translations into the amino acid realm, and therefore decrease the usefulness of the alignment. Plus, codon alignment can reveal synonymous vs nonsynonymous mutations. 

If you want to see an visualization of this application of alignment,[head over here](https://flea.ki.murrell.group/view/P018/sequences/). The sequences shown here are from an HIV envelope. If you click on on amino acid index, you can see a graph showing how the amino acid in that position evolved over time. You can explore the site - I recommend the tree section because it is pretty. 

I made a fake fasta of sequences that need to be codon aligned over at ```/srv/WI20/not_aligned.fasta```. Here is my step by step guide on making your own 

### 1. Align with mafft and replace initial gaps with N's

<details>
  <summary>Want to know the biological reason why we do this? click here</summary>
  
Sequencing starts at many different points, and we don't know ahead of time where the points are. If we don't identify where one sequence starts relative to another, we cannot begin comparing them. Why the N's? That's because the biological meaning of N's is different than that of gaps. Gaps in an alignment indicate deletions or insertions from one sequence to another, we predict that something was removed or added. In the case of different starting points, we know that *something* is supposed to be there, since we have information from other sequences. Thus it is not an insertion or deletion, but unkown nucleotides. We represent this with N's. 

</details></br>

  ⋅⋅**a.** Align the file with mafft. Syntax for this is provided in [lesson 8](https://github.com/ubicucsd/bootcamp_19_20/blob/master/8_Alignment1.md)
  
  <details>
  <summary>Can't figure out the correct syntax to call mafft? click here</summary>
  
```python
subprocess.call(["mafft", "--out", "nuc_aligned.fasta", in_file])
````
</details></br>

  ⋅⋅**b.** Go through each of the sequences in the mafft aligned file, see part 4 for how to do this
  
  ⋅⋅**c.** Count the number of initial gaps on each sequence
  
  ⋅⋅**d.** Degap each sequence and place it into a Seq object
  
<details>
  <summary>Running into type issues? This one is a bit of a pain, so let me give you this</summary>
  
```python
#This line converts the sequence to a string, replaces gaps with empty strings, and placed the result into a Seq object
sequence=Seq.Seq(str(seq_record.seq).replace("-", ""))
````
</details></br>
  
  ⋅⋅**d.** Prepend each sequence with n's. The number of n's should be equal to the number of initial gaps that sequence had.

#### **2.** Cut all of the sequences in such a way that they can be broken down into codons. 

<details>
  <summary>Usually this step is a bit more complicated. To learn more, click over here.</summary>
  
  
Even after accounting for the varying starting points in the sequence, we still have the issue of the reading frame. What if every single sequence starts at a nonsense location? All of the translated animo acids will be useless. A robust way of making the choice between starting each sequence from the first, second, or third nucleotide by seeing which one results in the longest total distance between stop codons in all of the sequences. In the interest of time, you are encouraged to skip this step and just find the which indices to pick in order to get mulitples of 3 in every sequence.

</details></br>


⋅⋅**a.** Go through each index from 0 to the number of sequences

⋅⋅**b.** Use the modulus to find whether the current sequences has a length divisible by three. What is a modulus? it's this ```%``` symbol. In math, this symbol will find the remainder for you. So 15%3=0, 16%3=1, 17%3=2, 18%3=0 and so on. 

⋅⋅**c** Based on the remainder you found in b, figure out how much of the sequence you should cut off the end. Syntax hint: ```seq[0:-1]``` will give you the sequence with the last nucleotide cut off. 


#### **3.** Translate the DNA into amino acids. This should be simple, I will leave this exercise up to you with one hint - google translating DNA with BioPython.

#### **4.** Mafft align the amino acid sequences from step 4. It's the same syntax as step 1, just put your animo acids in a file. 


#### **5.** You now have a good looking amino acid alignment, but how does it look as actual *codons*? We need to backtranslate from amino acids into codons by inserting the gaps from the amino acid alignment into the nucleotide codons themselves. 

⋅⋅**a.** Go through each amino acid in the aligned amino acid array

⋅⋅**b.** If it's a regular amino acid, go ahead and take next three nucelotides from the degapped nucleotide array (the one with the n's). If it's a gap, think about how many gaps that corresponds to in nuceotide space (hint: it's three) and insert those gaps. 







