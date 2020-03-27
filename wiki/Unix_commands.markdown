---
title:  "UNIX commands for computational chemistry "
date:   2019-07-03
layout: "git-wiki-post"
---

The most important UNIX commands needed for computational research

### About this tutorial

As a computational chemist, you should become acquainted with the basics of the Linux operating system. There is a wealth of free information about Linux available online and in several books. This tutorial is a starting point to help you begin using Linux and to understand the common commands useful to scientific computing.

### man (manual)

The easiest way to get more information on a particular Linux command or program is to use the *man* command followed by the item you want information on:

{% highlight git %}
$ man [program-name]
{% endhighlight %}

### ls (list)

{% highlight git %}
$ ls
{% endhighlight %}

list all files within a directory. Using the -l modifier gives more information such as read/write permissions. The -a option includes entries starting with "." The command can be combined with search strings in order to produce partial listings. For example

{% highlight git %}
$ ls test*.com
{% endhighlight %}

will list all files starting with the string "test" and ending with ".com".

### cd (change directory)

{% highlight git %}
$ cd [path/name]
{% endhighlight %}

change to the directory called "name"; when given without any arguments, the shell moves to the user's home directory. In order to change to a specific directory, the full pathname must be given e.g. cd /usr/local/bin changes to the directory /usr/local/bin. Moving up one directory can be done more easily by using

{% highlight git %}
$ cd ..
{% endhighlight %}

To display the full path of the current directory type

{% highlight git %}
$ pwd
{% endhighlight %}


### mkdir (make directory)

{% highlight git %}
$ mkdir [name]
{% endhighlight %}

creates directory "name" as a subdirectory of the current working directory.

### cp (copy a file)

{% highlight git %}
$ cp [file1] [file2]
{% endhighlight %}

copies "file1" to a new file called "file2"

### mv (moving a file)

{% highlight git %}
$ mv [file1] [file2]
{% endhighlight %}

changes the name of "file1" to  "file2"

### rm (deleting a file)

{% highlight git %}
$ rm [file]
{% endhighlight %}

removes file "file1". To remove a directory with all of its contents use the flag -r for recursive

{% highlight git %}
$ rm -r [directory]
{% endhighlight %}

Caution: deleting a file with *rm* cannot be undone. Make sure you have backups!

### find

 {% highlight git %}
 $ find -iname [file]
 {% endhighlight %}

finds file file in the current working directory as well as all of its subdirectories. The command

{% highlight git %}
$ find -iname "*400"
{% endhighlight %}

will thus locate all files whose names end with the string 400.

### grep

{% highlight git %}
$ grep string1 file2
{% endhighlight %}

search for string1 in file file2. For example, the command

{% highlight git %}
$ grep "B3LYP" test*.log
{% endhighlight %}

will search for string *B3LYP* in all files starting with test and ending with .log. Please observe that UNIX is case sensitive (e.g. "B3LYP" is different from "b3lyp"). If the string occurs repeatedly in one file, all occurrences are listed. Searching through all files of a given file system can be achieved with the -r modifier. The command

{% highlight git %}
$ grep -r /scr6/user99 -e "573442"
{% endhighlight %}

searches for number sequence 573442 in all files located in file system /scr6/user99. For searching binary files, use grep -a.

### top

{% highlight git %}
$ top
{% endhighlight %}

dynamic real-time view of the processes running on the system. Each task has a unique PID (process ID). To exit hit *control c*

### pwd (print working directory)

{% highlight git %}
$ pwd
{% endhighlight %}

prints the name of the current working directory.

### Input and Output (I/O)

 To save the output (stdout) from a program to a file use the redirection operator >

 {% highlight git %}
 $ example_program > my_output.txt
 {% endhighlight %}

 Input can  be given to a command from a file instead of typing it in the shell by using the redirection operator <

 {% highlight git %}
 $ my_command < programinput
 {% endhighlight %}

Alternatively, you can use the "pipe" operator:

{% highlight git %}
$ cat programinput | mycommand
{% endhighlight %}

Using the pipe operator, you can link commands together. The pipe will link stdout from one command to stdin of another command. In the above example we use the cat command to print the file to the screen (stdout), and then we redirect that printing to the command mycommand.

### ssh (secure shell)

SSH, or Secure Shell, is a remote protocol that allows users to access remote servers over the Internet. This allows you to remotely access a machine and execute shell commands in the same manner as you would if you were physically operating the remote computer.

{% highlight git %}
$ ssh {user}@{host}
{% endhighlight %}

The SSH key command instructs your system that you want to open an encrypted Secure Shell Connection. {user} represents the account you want to access. You will be prompted to enter the password for the requested account. Note: When you are typing your password, nothing will appear on the screen, but your password is, in fact being transmitted.

To display graphics (for example from plotting programs) you need to enables trusted X11 forwarding with the -Y option

{% highlight git %}
$ ssh -Y {user}@{host}
{% endhighlight %}
