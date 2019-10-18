# Hello World: What's the purpose of bioinformatics?

### Welcome to the Crash Course!

These activities are designed for students with zero to some programming experience and zero to some bioinformatics experience.  The crash course is always in development, so we ask for everyone to be patient and tell us when we are screwing up. Thanks!

## Big Bioinformatics Problems

As you progress through bioinformatics, you might find yourself developing technical skills without knowing what problems to tackle. In fact, there are several common and difficult problems in bioinformatics worth knowing. Therefore, this lesson will introduce you to of some of the things we, as bioinformaticians, do. 

*Note: In this course, we try to avoid lecture-based material. However, in this first lesson, we want to set a foundation of major problems bioinformaticians work on. We'll keep it short, and move on to the interactive part: working in a Linux environment!*

### I. The Alignment Problem ([full lesson](https://github.com/sabeelmansuri/Bioinformatics-Crash-Course/blob/master/5_Alignment.md))

The first problem we'll look at is alignment. **Alignment is the process of arranging sequences in a way to identify regions of similarity that may be a consequence of functional, structural, or evolutionary relationships between the sequences.** In order to understand why alignment is hard and practical problem, we will start by learning a bit about sequencing DNA. 

**Illumina Sequencing:** 

This is the most widely used "golden standard" method of sequencing DNA. The main idea is that DNA is fragmented into pieces, attached to a flow cell, copied many times over by PCR, and the complement strand is determined one nucleotide at a time. The modified nucleic acids Illumina uses during the generation of complementary strands emit light when they are bound. Illumina sequencing adds one base, measures the color output of the flow cell, adds a different base, measures the color output and repeats over and over. Here is a nice graphic depicting the process:

![graphic explaining illumina](http://www.3402bioinformaticsgroup.com/wp-content/uploads/2016/07/NGS.png)

There is also a video explaining Illumina sequencing which you will watch in *every single Bioinformatics class* you take in the future, found [here](https://www.youtube.com/watch?v=fCd6B5HRaZ8). 

Want to learn more about different sequencing methods? We have a document describing non Illumina methods in more detail [over here](https://github.com/sabeelmansuri/Bioinformatics-Crash-Course/blob/master/extra_sequencing.md)

***Comparing Sequencing Methods and Why they Matter***

What do all sequencing method have in common(yes, even PacBio)?  All produce ridiculous quantities of small DNA sequences. Illumina produces  300 million to 4 billion reads per run, with a selection of read lengths ranging from 50 base pairs to 300 base pairs. Meanwhile, Sanger produces 50000 sequences at lengths varying from 800 to 1000 base pairs. To give some perspective, the typical animal of interest is a human and those have 3.0×10^9 base pairs. Individual human genes range from 1148 to 37.7 kb (average length = 8446 bp,s.d. = 7124). 

These tiny reads overlap all over the place. If you imagine the true sequence these reads came from and place the reads where they came from, you will get many reads piled up over every base pair in the true sequence. The more reads pile up, the more accurately you can predict the actual sequence. A common measure that rates the robustness of an alignment is coverage:

![coverage](https://slideplayer.com/slide/5083621/16/images/4/Definition+of+Coverage.jpg)

Besides helping us reconstruct the DNA put into a sequencer, alignment can also
- Determine which organism an unkown sequence comes from
- Pinpoint locations where mutations have occured relative to a reference sequence
- Help determine the evolutionary distance between sequences

TLDR: There are numerous complex applications of bioinformatics algorithms, from functional structure predictions to ancestral reconstructions. Alignment serves as the foundation for many of these algorithms, making basic sense of the incomprehensible mass of DNA that sequencing gives us. 

### II. The Clustering Problem ([full lesson](https://github.com/sabeelmansuri/Bioinformatics-Crash-Course/blob/master/6.%20Clustering.md))

**Clustering is the the process of assigning data points to groups in such a way that the elements in a group/cluster are more similar to each other than they are to those in other groups.** The definition of what it means to be similar can vary and is determined by the function we use to measure distance between two points. 

#### Some types of clustering:
- Hierarchical Clustering: Repeatedly combines the closest points into a cluster that is the hybrid location of both points. The reason this is interesting to us is that it forms a tree of clusters, which can represent a tree of related genes which can be used to infer homology.

![hierarchical cluster](https://upload.wikimedia.org/wikipedia/commons/a/ad/Hierarchical_clustering_simple_diagram.svg)

- K-means Clustering: Separates a dataset into k groups of points in such a way that the members of a cluster are as close as possible to the center of the cluster they belong to. This type of clustering can check that the datapoints we are observing cluster together by tissue type, experimental conditions, time points, etc. 

- Fuzzy Clustering: Datapoints are not definitively assigned to a specific cluster, rather they are given a likelihood of belonging to a cluster. This can be used to ascertain levels of co expression between genes, revealing genes which may be under common regulatory control. 

## The Linux Terminal

## Getting Set Up for ssh

Bioinformatics is often space and computation intensive, so we outsource our computational work to a bigger computer called a server.

***Why Linux?*** The majority of servers run on Linux, a free operating system which inherited its predecessor's (GNU's) mission to give users freedom. Linux is completely open source, allowing users to see and modify any part of its inner workings. Linux is also extremely stable, allowing servers to be up for years at a time without restarting. 

## Task 1: get an ssh client 
***Secure Shell(ssh):*** a protocol which creates a secure channel for two computers to communicate even over an unsecured network. This is how we will connect to EC2. 

**How to prepare for ssh:**

**Windows:** I do not like Windows terminals, you do not like Windows terminals, no one likes Windows terminals. Go to [The Putty website](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html) and download Putty. A GUI should appear to guide you through the installation. 

**Mac or Linux:** No preparation necessary, since you already have a native ssh client. 

## Task 2: Aliview

Okay, we have had enough of conceptual stuff. Let's get at it with some cool visuals. 

You will need java to proceed. If you do not have it, go [here](java.com) and press the big red install button in the middle of the page. 

Aliview is a sequence viewer with a bunch of built-in tools, including alignment tools. We will use Aliview to see what a typical dataset looks like coming out of the illumina sequencing machine and what it means, visually, to align the sequences. Click [here](http://www.ormbunkar.se/aliview/) and go to download the stable version for your OS. Next, download a neat dataset I have for you from [here](https://drive.google.com/a/ucsd.edu/file/d/1iIzDwKm2k_VOnen0BRwtZ9Jlf81h3t5g/view?usp=sharing). Launch aliview, click file->open file->PC64_V48_small_reconstructed_seqs.fasta. Scroll to the right and notice the mess that begins to form as you scroll. These sequences are sourced from the same gene and have gone through many steps, so the differences between them are quite likely to be real. Click align in the upper left corner and click realign everything. Now scroll forward and observe the gaps that have been inserted by the aligner.

These sequences are actually from HIV-1 glycoprotein envelopes from a person with broadly neutralizing antibodies against HIV-1(sequenced with PacBio). An effective HIV-1 vaccine should evoke the production of broadly neutralizing antibodies, so it is important to study the structure of the envelope that caused these antibodies to develop in individuals. The steps leading up to the data you downloaded allowed us to reveal a few specific strains. Now that we have them aligned, we can start asking questions about the differences between them and their evolution(which is important to figuring out how to counteract them). If you want to go into the nitty-gritty biology behind these sequences, go [here](https://www.cell.com/immunity/pdf/S1074-7613(17)30479-X.pdf). I will make a reconstruction of the evolutionary relationship for you to look at. 

The SNP differences are obvious, but you will notice that there are weirder differences - like the >100 bp gaps formed in the middle. The truth is that I do not know why those are there! Those can be 

1. Real differences in hypervariable loops
2) non-functional virions
3) RT/PCR artifacts

Here's another cool tool: Blast will quickly look up a sequence in the NCBI database and spit out similar sequences it finds. Now, go to aliview and click edit->delete all gaps in all sequences. Copy the first sequence into the clipboard. Google NCBI Blast and open up the first result. Paste the sequence into the big box in the top of the page and click **BLAST** at the bottom. After a few seconds, Blast should link you to a whole lot of glycoprotein that are... from the article these were published in! This is one of the basic uses of Blast - to figure out where a sequence comes from. 

## Task 3: Explore your EC2

In your email, you should have a password from me. Your username for EC2 is the same as your ucsd username.

**Windows:** Open putty, paste ```ec2-13-59-255-161.us-east-2.compute.amazonaws.com``` into the Host Name section and select port 22 and SSH on that same page. Type a name under Saved Sessions and click the Save icon on the right. Now, press Open and type your username and password when prompted. 

**Mac or Linux:** Right click anywhere and click open terminal. You should see a prompt that looks something like  ```mchernys@mchernys-ThinkPad-T430:~/Desktop$```. Next, copy paste this command into the terminal and press enter ```ssh your-username@ec2-13-59-255-161.us-east-2.compute.amazonaws.com```. Note: use Ctrl-shift-V to paste into terminal. Please replace "your-username" with your actual username. Save the command you used somewhere so you can copy paste it in the future. 

You probably have your own password in mind for your account. Type ```passwd username``` and follow the prompts to set your own password. 

Discover your identity. Type `whoami` into the window that just opened up and hit `enter`. And just like that you're talking
with your computer, you bioinformatician, you.

## How do I see what a command does?

Anytime you need a refresh on what a command does, type the command line with the --help option like so: ```ls --help```. If that does not work, try ```man ls```. I will go over why different commands have different help syntax in a bit. 

## Navigation, manipulation, and permission

In order to get started, we need to be able to do the same thing we do in a file explorer in the command line. You may find it inconvenient at first, but with time these commands become faster and more versatile than the file explorer's interface. 

The forward slashes in a terminal console represent directories, with the home directory being a ```~```. Your default folder on EC2 is your use folder, which is ```~/username```. This means the folder named after your username is a subfolder of the home folder, which is represented by ```~```. 

```cd```(change directory) Type cd followed by the directory's path to navigate a terminal to that directory. ```.``` is current directory and ```..``` is the parent of the current directory. 

```ls```(list files) prints out the contents of a directory. There are tons of options for this command - my favorite is ```ls -lah``` , since it prints the directory contents in list format(```-l```), includes hidden files/folders(```-a```), and makes the storage sizes more readable for humans(```-h```). 


```mkdir```(make directory) Creates a directory with the same name as the argument you give it. 

---

### TODO: Make a Software Folder

Navigate your terminal to your home directory(the directory named after your UCSD username) using ```cd```. Type ```mkdir software``` and press enter. Type ```ls``` to see the changes you have made. The reason for a software folder is to keep your software in it, oddly enough. Usually, you would place executables in the /bin system folder, but you are not the admin so you cannot access that folder :( . This is often the case when you ssh into a system, so get used to having a dedicated software folder.  

---

3 tasks:

1. Print all of the contents of the downloaded dataset to terminal window ("standard output")
2. Print how many lines there are in the file
3. Print how many lines there are in the file **THAT CONTAIN GENETIC DATA** (no headers)
