# Advanced Command Line + Python

### Welcome Back!

Welcome to week 2 of the Crash Course. Today's lesson will be split into two distinct topics. The first builds on the skills you developed [last session](/1_Welcome.md), testing and enforcing your understanding of the command line. The second will introduce you to Python.

## Jumping in the Deep End

Last session, we said that we prioritize interactive lessons. We're going to waste no time living by this mantra.

## Getting yer feet wet

Here are a quick batch of tasks/exercises using the command line (and different commands) in roughly increasing order of difficulty:

##### (0.) Copy-paste the following two commands into your terminal in order and hit enter.
`cd ~`  
`cp -r ../smansuri/week2/ ~`

##### 1. Enter the directory that was just created

##### 2. Name all of the files (not other directories!) inside this directory. How many are there? (hint: use `ls`)

##### 3. How many lines are in file1.txt? (hint: use `wc`)

##### 5. How many lines inside file2.txt contain the characters (not the word) "no"? This means "anNOunce" will count too. (hint: use a pipe of `grep` and `wc`)

##### 6. (Challenge 1) Execute the "instructions" file. Follow the instructions. (hint: `./instructions`)

##### 7. (Challenge 2) Attempt to delete the directory you downloaded in step 0. This directory is called "parent".


## Downloading

```scp```(secure copy) is a command used to copy files from one machine to another. The first argument is the source location, while the second argument is the destination. ```scp file.txt my_username@dns_address.com:/home/my_username/docs```

```curl``` Will download stuff for you. The most simple and relevant combination of options is ```curl -L https://examplelink.com -outdir .``` which will download from https://examplelink.com into the current directory (indicated by the dot). 

```apt-get```Handles packages from the apt library for Debian based systems. However, this installs packages system-wide so you are not going to be able to use it on EC2. The mac equivalent is homebrew. ```sudo apt-get install google-chrome-stable``` will install chrome. 

---

### TODO: Get FastQC and Kallisto

Quality of genetic information is important! FastQC is the gold standard for quality control in the bioinformatics field. Google "download FastQC" and find the instructions. Make sure to select the correct package for a linux system, then go think about how you would get that package onto EC2. There are two main ways to do this. 

*Hint: the last three commands mentioned contain both of the two ways you can get FastQC onto EC2*

Next, we need Kallisto, which describes itself as "a program for quantifying abundances of transcripts from RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads." Google "download kallisto" and find the appropriate file (it should be a .tar.gz). Use the same method you used for FastQC to transfer it to EC2(or challenge yourself to find the second way).

---

## Unpackaging

Much of the data people want to download is large, but they want it fast. That's why things like .zip, .tar, .gz and such exist. Those are the file extensions of compressed data. In order to make software work, it must be unpackaged.

```unzip``` Is exactly what it sounds like. This command unzips .zip file types. 

```tar```(tape archive) Is the command linux uses to package and unpackage stuff. This command has an incomprehensible amount of confusing options, so let me just copy paste the ones you should care about. ```tar -xvf file.tar.gz -C .``` unpacks a .tar.gz file into the current directory and ```tar -xvf -C .``` unpacks a .tar file into the current directory. The -C option indicates the files' destination.

---

### TODO: Unpackage your FastQC Kallisto

Now you have your fastqc*.zip in your software folder. In order to use it, you're going to have to unpackage it. Look a couple lines up to figure out how. 

The same goes for your kallisto .tar.gz. This is a filetype you will run into often when dealing with linux, since it is the default way linux compresses a folder. If you look in the Unpackaging section of this document, you will find instructions on how to open up this strange creature. 

---

## Compilation

What is compilation? It is the conversion of one programming language into another. Typically, it is a conversion of what's known as a high-level language (C, Java, Python, etc) to a low level language (binary, assembly). CPUs understand only very very very basic logic, so a super smart program called a compiler has to convert your convoluted and messy code into the simple delicious porridge that the CPU can eat(execute).

**shell scripts and python files** do not need compilation.

**java** compiles by ```javac filename.java```

**C** ```gcc -c filename.c``` to compile and assemble. 

Note: Much of the time, software you download online is already in binary form so there is no need to compile. This is not always the case!

## Execution

**/bin**(binaries) contains your executable files and shells. The computer has a list of folders it searches through to find executables when you type a command and the /bin directory is one of them. When you download software, you should place the executable file or a symbolic link into the /bin directory.

**shell script** a simple ```./executable``` will suffice to execute a script. 

**python** will automatically compile for you before executing with the command ```python filename.py```. 

**java** can be executed with ```java compiledfilename```

**C** is executed like a shell script ```./out```

---

### TODO Actually Quantify Stuff

Okay, we now how to execute now. The FastQC folder you have now contains an executable called fastqc. Go into the folder containing the executable, type ```./fastqc --help``` to see usage instructions. Your task is simply to run the fastqc on each of the files sitting in the ```/srv/lesson2/sub_Gastric_Ctr``` and ```/srv/lesson2/sub_Gastric_Affect```directories. FastQC will produce some .html files, which I will just show on the large screen and explain in the interest of saving time. 

Kallisto is a little more complicated. Our first step is to build a kallisto index file, which will assign an index to each RNA transcript we will quantify and optimizes the quantifying procedure in general. Where did we get these RNA transcripts you ask? Good question! I googled "mouse lncRNA database" and eventually happened upon a [website](https://www.gencodegenes.org/) where there was a fantastically easy to download .fasta file full of lnc RNA sequences from mice. Anyways, go back to the kallisto website and look at the instructions on how to index a file. The list of target sequences I got from gencode is at ```/srv/lesson2/gencode.vM17.lncRNA_transcripts.fasta``` and you might want to specify the name for the destination file. 

Next, you will need to look at the ```kallisto quant``` command. Specify an output folder, the index you made just now, set --bootstrap-samples=100, and finally include the forward and reverse paired end files at the end. You will run the quant command 4 times for the 2 control and 2 affected files. You might want to create 4 separate output folders in order to keep everything organized. 

The last part will be a bit of a walkthrough, since it is kind of complicated (I don't understand every option either, don't worry about it). Open up R  by simply typing R into the terminal and pressing enter. Below is the template for what you need to execute in R in order to compare transcription levels. Note the parts where you need to enter a path and replace them with your own filepaths


```
library("sleuth")
#paths to Kallisto outputs from MPNST study
sample_id = c("/home/mchernys/Documents/UCSD_Classes/Spring_2018/CSE_185/final_project/kallisto_output/Gastric_Ctr_Rep1", "/home/mchernys/Documents/UCSD_Classes/Spring_2018/CSE_185/final_project/kallisto_output/Gastric_Ctr_Rep2", "/home/mchernys/Documents/UCSD_Classes/Spring_2018/CSE_185/final_project/kallisto_output/Gastric_Affect_Rep1", "/home/mchernys/Documents/UCSD_Classes/Spring_2018/CSE_185/final_project/kallisto_output/Gastric_Affect_Rep2")

kallisto_dirs = file.path(sample_id)
s2c = read.table(file.path("/home/mchernys/Documents/UCSD_Classes/Spring_2018/CSE_185/final_project/sleuth_Gastric_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c = dplyr::mutate(s2c, path = kallisto_dirs)

so = sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)


write.table(sleuth_significant, "/home/my_username/sleuth_output/", sep="\t", quote=FALSE)
gg=plot_transcript_heatmap(so, transcripts, units = "tpm", trans = "log", offset = 1)
ggsave("/home/my_username/plots/", plot=gg)
```
