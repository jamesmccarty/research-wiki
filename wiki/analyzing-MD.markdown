---
title:  "Analyzing a MD Trajectory"
date:   2019-07-03
layout: "git-wiki-post"
---

Analyzing a MD Trajectory

### About this tutorial

This tutorial deals with extracting basic information from a MD trajectory performed in GROMACS. If you have not yet completed the [Basics of a GROMACS simulation](2019-08-01-gromacs), you should do that tutorial before proceeding to this tutorial. This tutorial will introduce you to a few basic tools to analyze a generated trajectory. For your own system, you should formulate for yourself some ideas about the types of data you will want to collect.

### Converting trajectory formats

The Gromacs command *trjconv* is a post-processing tool to manipulate the atomic coordinates, extract a subset of the coordinates, correct for periodic boundary conditions (pbc), or manually alter the trajectory (time units, frame frequency, etc.). Here we will correct for the periodicity of the system. The protein will diffuse through the unit cell and may appear to "jump" across to the other side of the box. To account for such actions, apply the following command:

{% highlight git %}
$ gmx trjconv -s md_trajectory.tpr -f md_trajectory.xtc -o md_trajectory_noPBC.xtc -pbc whole -ur compact
{% endhighlight %}

where -s signals the input .tpr file that was used to run the MD simulation (here called *md_trajectory.tpr*) and -f signals the trajectory that was produced by the GROMACS *mdrun* command. Here we are using the compressed .xtc format. The -o is the output modified trajectory file. The -pbc whole flag signals to make molecules that were split by the periodic boundary conditions whole and the -ur compact flag signals to keep all molecules in the original unit cell. When prompted select option 0 for the entire system.

In the remainder of this tutorial we will work with this "corrected" trajectory *md_trajectory_noPBC.xtc*

If you want to look at a movie of the whole trajectory, you can convert the xtc file to a pdb file using *trjconv* command

{% highlight git %}
$ gmx trjconv -s md_trajectory.tpr -f md_trajectory_noPBC.xtc -o md.pdb
{% endhighlight %}

This will convert the trajectory into .pdb format which can be viewed in PyMOL or VMD.

### Root-mean-square deviation (RMSD) of atomic positions

The root-mean-square deviation of atomic positions (or simply root-mean-square deviation, RMSD) is the measure of the average distance between the atoms (usually the backbone atoms) of superimposed proteins. This will give a measure of the overall change in the conformation of the protein as the rmsd with respect to a reference state. We can use the first frame of the trajectory as a reference state. First we will extract this frame using *trjconv*

{% highlight git %}
$ gmx trjconv -s md_trajectory.tpr -f md_trajectory.xtc -o frame1.pdb -pbc whole -ur compact -dump 0
{% endhighlight %}

where we again correct for periodic boundary conditions as before, but here we output just the first frame with the flag -dump 0.

To compute the RMSD enter the command:

{% highlight git %}
$ gmx rms -s frame1.pdb -f md_trajectory_noPBC.xtc -o rmsd-vs-start.xvg
{% endhighlight %}

First the algorithm will align each structure to the reference configuration (frame1.pdb). We must select which atoms to use for the least-squares aliment. Common choices include the backbone atoms or all heavy atoms. When prompted, select 4 for Backbone atoms. The next prompt will ask us to select a group for the RMSD calculation. You can again experiment with different selections. Here we will select 4 again for the Backbone atoms.

Plot the generated rmsd-vs-start.xvg using the Xmgrace plotting tool

{% highlight git %}
xmgrace rmsd-vs-start.xvg
{% endhighlight %}

A plot of the RMSD is shown below:

![]({{ site.url }}{{ site.baseurl }}/images/rmsd-backbone-vs-start.png){: style="width: 600px; border: 10px"}

An alternative metric is to compute the root mean square fluctuatoins (RMSF). The RMSF captures, for each atom, the fluctuation about its average position. This gives insight into the flexibility of regions of the peptide. The RMSF (and the average structure) are calculated with the *rmsf* command. We are most interested in the fluctuations on a per residue basis, which is controlled by the flag -res.

{% highlight git %}
$ gmx rmsf -s frame1.pdb -f md_trajectory_noPBC.xtc -o rmsf-per-residue.xvg -ox average.pdb -res
{% endhighlight %}

Have a look at the graph of the RMSF with xmgrace and identify the flexible and rigid regions of the protein.

{% highlight git %}
xmgrace rmsf-per-residue.xvg
{% endhighlight %}

A plot of the RMSF for ubiquitin is shown below:

![]({{ site.url }}{{ site.baseurl }}/images/rmsf-per-residue.png){: style="width: 600px; border: 10px"}

A side product of the RMSF calculation is the average protein structure over the course of the simulation. Note that the average protein structure is not necessarily a physically relevant structure if there are large conformational changes during the simulation. To get a better measure of the convergence of the simulation to equilibrium, we compute the RMSD again, but this time using the average structure as the reference structure.

{% highlight git %}
$ gmx rms -s average.pdb -f md_trajectory_noPBC.xtc -o rmsd-vs-average.xvg
{% endhighlight %}

Again select 4 for the Backbone atoms for both the alignment and RMSD calculation.

Plot the generated rmsd-vs-average.xvg using the Xmgrace plotting tool

{% highlight git %}
xmgrace rmsd-vs-average.xvg
{% endhighlight %}

A plot of the RMSD is shown below:

![]({{ site.url }}{{ site.baseurl }}/images/rmsd-backbone-vs-average.png){: style="width: 600px; border: 10px"}

At equilibrium about a single stable configuration, the RMSD value should be level and fluctuate about the mean.

### Radius of gyration

The module gyrate allows you to check the radius of gyration of the system. If the protein is unfolding or adopting "open" configurations, the radius of gyration will increase.

{% highlight git %}
gmx gyrate -s frame1.pdb -f md_trajectory_noPBC.xtc -o gyrate.xvg
{% endhighlight %}

Plot the generated gyrate.xvg using the Xmgrace plotting tool

{% highlight git %}
xmgrace gyrate.xvg
{% endhighlight %}

A plot of the radius of gyration is shown below:

![]({{ site.url }}{{ site.baseurl }}/images/gyrate.png){: style="width: 600px; border: 10px"}
