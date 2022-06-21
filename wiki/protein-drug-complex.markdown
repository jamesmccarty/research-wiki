---
title:  "Protein-Ligand Complex"
date:   2019-07-03
layout: "git-wiki-post"
---

## Protein-Ligand Complex ##

CDK5/(R)-roscovitine Complex

### About this tutorial ###

In this tutorial we will use the Antechamber package of [AmberTools](https://ambermd.org/AmberTools.php) to create topology and structure files that can be used in a [GROMACS molecular dynamics simulation](2019-08-01-gromacs). The [Generalized Amber Force Fields (GAFF)](http://ambermd.org/antechamber/gaff.html) are a general set of force field parameters to describe small organic molecules such as drugs. In this tutorial we will parameterize the small molecule (R)-roscovitine in complex with the cyclin-dependent kinase, CDK5. The starting pdb file can be downloaded from the Protein Data Bank (PDB code: [1UNL](https://www.rcsb.org/structure/1UNL)).

Software needed to complete this tutorial:

[GROMACS](2019-08-01-gromacs) <br/>
[AmberTools](https://ambermd.org/AmberTools.php) <br/>

### Prepare starting pdb files ###

The pdb file 1UNL contains the p25 activator, the CDK5 protein, and the drug (R)-roscovitine. We need to prepare two separate topologies: (1) the protein CDK5 topology using GROMACS *pdb2gmx* and (2) the ligand topology using AmberTools.

First, we need to extract just the chain A coordinates:

{% highlight git %}
grep ATOM 1unl.pdb | grep ' A ' > cdk5.pdb
{% endhighlight %}

And we will also extract the (R)-roscovitine coordinates:

{% highlight git %}
grep RRC 1unl.pdb > roscovitine.pdb
{% endhighlight %}

### Prepare protein topology file ###

To generate the protein topology file we use the GROMACS *pdb2gmx* tool

{% highlight git %}
gmx pdb2gmx -f cdk5.pdb -o cdk5.gro -ignh -p cdk5.top
{% endhighlight %}

When prompted, choose (6) AMBE99SB-ILDN for the force field, and (1) TIP3P for the water model.
The file *cdk5.gro* contains the coordinates (with hydrogens added) of the CDK5 protein, and the file *cdk5.top* contains the topology and force field parameters for the protein and water molecules. Now we must add the ligand, but the parameters for the ligand are not included in the protein force fields. We must add these parameters ourselves. In this tutorial, we will use the Generalized Amber Force Fields (GAFF) for small molecules.

### Prepare a parameter and coordinate file for (R)-roscovitine ###

Make sure your AMBER_PREFIX is set by typing in the termimial:

{% highlight git %}
$ source /path/to/amber18/amber.sh
{% endhighlight %}

where /path/to/amber18 is the path to ambertools on your machine. First, we use the *reduce* command to add hydrogens to the pdb file for (R)-roscovitine.

{% highlight git %}
$ reduce roscovitine.pdb > roscovitine_h.pdb
{% endhighlight %}

Our next step is to use antechamber to create a *mol2* file which is required to build our topology file.

{% highlight git %}
$ antechamber -i roscovitine_h.pdb -fi pdb -o roscovitine_h.mol2 -fo mol2 -c bcc -s 2
{% endhighlight %}

Here the -i signals the input file and -fi specifies that the input file is in pdb format. The -o specifies the name of the output file to write and we use -fo to specify that we want the output to be written in mol2 format. The -c bcc option tells antechamber to use the AM1-BCC semi-empirical charge model to calculate the partial charges on each atom. Normally, we would calculate partial charges from the output of a quantum chemistry calculation using the restrained electrostatic potential (RESP) procedure, but for the purpose of this tutorial the semi-empirical AM1 model is sufficient. The -s 2 defines the verbosity - how much information antechamber will print to the screen during the procedure.

A number of new files are produced by this program. The *sqm.xxx* files are the input and out files for the sqm quantum chemistry code used by Antechamber to compute the point charges on each atom. Check to make sure the procedure completed with no errors by looking at the end of the sqm.out file. The last line should read

{% highlight git %}
---------- Calculation Completed ----------
{% endhighlight %}

The output mol2 file, *roscovitine_h.mol2*, contains the definition of our drug including all of the charges and atom types. The GAFF parameters are located in $AMBERHOME/dat/leap/parm/gaff.dat. At this stage we want to check to make sure there are no missing parameters in our molecule. We can use the utility *parmchk* to test if all the parameters required are available

{% highlight git %}
$ parmchk2 -i roscovitine_h.mol2 -f mol2 -o roscovitine_h.frcmod
{% endhighlight %}

The above command produces the output *roscovitine_h.frcmod* file which is a parameter file that we can later load to add missing parameters. It will contain all of the missing parameters. Antechamber will try and fill in these missing parameters by analogy to similar structures. Looking at our frcmod file we see that 8 improper angles had missing parameters. For the sake of this tutorial, we will assume that the parameters Antechamber has suggested for us are acceptable. If you see any parameters listed with the comment "ATTN: NEEDS REVISION" then it means that Antechamber could not determine suitable parameters and so you must manually provide these before you can proceed with the simulation.

Next we load the amber software LEaP and load in the GAFF force field parameters.

{% highlight git %}
$ tleap -f oldff/leaprc.ff99SB
{% endhighlight %}

followed by
{% highlight git %}
$ source leaprc.gaff
{% endhighlight %}

Now we can load our drug molecule (roscovitine_h.mol2)

{% highlight git %}
$ roscovitine = loadmol2 roscovitine_h.mol2
{% endhighlight %}

We now need to load our frcmod file in order to tell tleap the how to treat the missing parameters.

{% highlight git %}
$ loadamberparams roscovitine_h.frcmod  
{% endhighlight %}

We can save a library file for the (R)-roscovitine molecule (roscovitine.lib) and a topology file prmtop and inpcrd file (roscovitine.prmtop, roscovitine.inpcrd)

{% highlight git %}
$ saveoff roscovitine roscovitine.lib
{% endhighlight %}

{% highlight git %}
$ saveamberparm roscovitine roscovitine.prmtop roscovitine.inpcrd
{% endhighlight %}

Now we can exit tleap
{% highlight git %}
$ quit
{% endhighlight %}

We now have the Amber files *roscovitine.prmtop* and *roscovitine.inpcrd* which contain the force field parameters and input coordinates in Amber format. We now need to convert these into a GROMACS topology and structure file using the Python interface ACPYPE

{% highlight git %}
acpype -p roscovitine.prmtop -x roscovitine.inpcrd
{% endhighlight git %}

You should now have the Gromacs files *RRC_GMX.gro* and *RRC_GMX.top* which contain the GROMACS topology and coordinates for (R)-roscovitine

### Combining Ligand and Protein Files ###

Now we have our protein files *cdk5.gro* and *cdk5.top* and the ligand files *RRC_GMX.gro* and *RRC_GMX.top*. In order to have one set of GROMACS files we can combine them with the python script [ligand2gro.py](https://github.com/jamesmccarty/Protein-Ligand-Complex/blob/master/ligand2gro.py).

{% highlight git %}
python ligand2gro.py -protein cdk5.gro -ligand RRC_GMX.gro -proteintop cdk5.top -ligandtop RRC_GMX.top
{% endhighlight git %}

The output files will be a GROMACS structure file for both the ligand and protein called *complex.gro* and a topology file *complex.top* containing the force field parameters.

You can now create a simulations box and add solvent molecules

{% highlight git %}
gxm editconf –f complex.gro –o complex_newbox.gro –c –d 1.0 –bt cubic

gmx solvate –cp complex_newbox.gro –cs spc216.gro –o complex_solv.gro –p complex.top
{% endhighlight git %}

We also add ions as before using the [ions.mdp](https://github.com/jamesmccarty/Protein-Ligand-Complex/blob/master/ions.mdp) file:

{% highlight git %}
gmx grompp –f ions.mdp –c complex_solv.gro –p complex.top -o ions.tpr -maxwarn 1

gmx genion -s ions.tpr -o complex_solv_ions.gro -p complex.top -conc 0.15 -neutral
{% endhighlight git %}

Finally, we perform an energy minimization with [em.mdp](https://github.com/jamesmccarty/Protein-Ligand-Complex/blob/master/em.mdp):

{% highlight git %}
gmx grompp -f em.mdp -c complex_solv_ions.gro -p complex.top -o em.tpr

gmx mdrun -v -deffnm em
{% endhighlight git %}

At the end of this run we now have an energy minimized structure of the protein/ligand complex called *em.gro*. From this structure file and with our topology file *complex.top* we can continue performing an NVT equilibration simulation as described in [Basics of a GROMACS simulation](2019-08-01-gromacs).
