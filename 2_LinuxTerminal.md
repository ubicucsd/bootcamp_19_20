# The Linux Terminal

Bioinformatics is often memory and computation intensive, so we'll be outsourcing our computational work to a bigger computer called a server. This server will run on the Linux operating system. (This is no different than your laptop running on a Windows or Mac operating system.)

***Why Linux?*** The majority of servers run on Linux, a free operating system which inherited its predecessor's (GNU's) mission to give users freedom. Linux is completely open source, allowing users to see and modify any part of its inner workings. Linux is also extremely stable, allowing servers to be up for years at a time without restarting.

## How Can I Connect?

***Secure Shell(ssh)*** is a protocol which creates a secure channel for two computers to communicate. This is how we will connect to our server.

**Mac or Linux users:** skip to **Exploring the Server**. Your device has a built-in ssh client!

**Windows users:** We do not like Windows terminals, you do not like Windows terminals, no one likes Windows terminals. Go to [The Putty website](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html) and download Putty. A window should appear to guide you through the installation.

## Exploring the Server

We've created an account for you on our server (which is, by the way, called EC2). Your username for EC2 is the same as your UCSD username. Your password has been emailed to you.

### Connecting

**Mac or Linux:** Open an application on your device called "Terminal". This should open up a window that prompts you to enter some text. Copy-paste this command into the terminal: ```ssh your-username@ec2-3-16-81-18.us-east-2.compute.amazonaws.com```. Replace "your-username" with your UCSD username, and press enter. You're now connected to EC2!

*Note: you may need to use Ctrl-shift-V to paste into terminal.*

**Windows:** Open putty, paste ```ec2-3-16-81-18.us-east-2.compute.amazonaws.com``` into the Host Name section and select port 22 and SSH on that same page. Type any name of your choice under Saved Sessions and click the Save icon on the right. Now, press Open and type your username and password when prompted. 

### Your First Commands

1. Discover your identity. Type `whoami` into the window that just opened up and hit `enter`. And just like that you're talking
with your computer, you bioinformatician, you.

2. The password we provided you with isn't particularly safe. Let's change your password using the terminal. Type ```passwd [your-username]``` and follow the prompts to set your own password. (Replace "your-username" with your username.)

### So... Commands?

Commands are things you can type into the terminal to perform different actions. There are an endless number of commands, each with a ridiculous amount of options, so **do NOT attempt to memorize them on the first try**. 

How do you know what command to use to do sometihng you want? Simple: for now, we'll explain commands as you need them. Actually learning (or memorizing) the commands will come naturally from repeated use of the terminal.

*Note: Anytime you need a refresher on what a command does, type the command line with the --help option like so: ```ls --help```. If that does not work, try ```man ls```. One (or both) of these will pull up information on how to use the command. Can you figure out what the ls command does?*

### What Am I looking At? What Can I Do? How Can I Do It?!

When you open up your laptop, you are presented with your "Desktop". When you open up a terminal (or connect to a server via ssh!), you are presented with your "Home".

When you want to open up the "Pictures" folder on your laptop, you find the folder labeled "Pictures", and then open it. Uh-oh... we don't know how to (1) find something or (2) open something!

Let's take a step back and talk about all the files in your computer are organized. As you know, you can have different files stored in different places on your computer. You do this by creating folders (inside of other folders), and creating files inside of them. This is **exactly how our server works**, except we call folders "directories".

Similar to how "Desktop" is a folder on your computer, "Home" is a directory on your account on the server. And similar to how your "Desktop" can have folders created in it, so can "Home". You can look at what's inside with the `ls` command, which is short for "list". Type `ls` to see what is in your "Home" directory. (We now know how to find something!)

The answer: nothing! Let's change that by creating a directory. You can do this with the `mkdir` command. Let's create a directory called "software" by doing `mkdir software`. Confirm that this worked by looking at your home directory again.

--

The forward slashes in a terminal console represent directories, with the home directory being a ```~```. Your default folder on EC2 is your use folder, which is ```~/username```. This means the folder named after your username is a subfolder of the home folder, which is represented by ```~```. 
  

### [Your Turn!] Investigate Genetic Data

3 tasks:

1. Print all of the contents of the downloaded dataset to terminal window ("standard output")
2. Print how many lines there are in the file
3. Print how many lines there are in the file **THAT CONTAIN GENETIC DATA** (no headers)


### Once you're done, show your answers to an instructor to get checked off. Congratulations! You've completed the first lesson of the Bioinformatics Crash Course!

## Commands Cheat Sheet

```ls```(list files) prints out the contents of a directory. There are tons of options for this command - my favorite is ```ls -lah``` , since it prints the directory contents in list format(```-l```), includes hidden files/folders(```-a```), and makes the storage sizes more readable for humans(```-h```). 

```mkdir```(make directory) Creates a directory with the same name as the argument you give it. 

```cd```(change directory) Type cd followed by the directory's path to navigate a terminal to that directory. ```.``` is current directory and ```..``` is the parent of the current directory. 




