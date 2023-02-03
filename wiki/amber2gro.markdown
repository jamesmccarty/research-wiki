---
title:  "Converting from Amber to GROMACS format"
date:   2023-02-03
layout: "git-wiki-post"
---

## Converting from amber format into GROMACS ##

### About this tutorial

Ambertools is a freely distributed program for preparing molecular dynamics simulations. If you build your system using Ambertools, then you can convert the Amber topology and coordinate files to GROMACS format using the ACPYE code [https://github.com/alanwilter/acpype](https://github.com/alanwilter/acpype).

This tutorial will go through the steps of converting an Amber topology and coordinate file into GROMACS format. In this example, we assume you have two files generated from Ambertools:
{% highlight git %}
system_solv.prmtop  -  Amber topology
system_solv.inpcrd  -  Amber coordinate file
{% endhighlight %}  

The necessary files to run this tutorial can be found [here](https://github.com/jamesmccarty/amber2gromacs).

### Convert amber topology and structure file to GROMACS topology and structure file using ACPYPE

Make sure acpype is in your path and type:
{% highlight git %}
acpype -p system_solv.prmtop -x system_solv.inpcrd -b system_solv
{% endhighlight %}  

Running this command will produce two files:
{% highlight git %}
system_solv_GMX.top  -  GROMACS topology
system_solv_GMX.gro  -  GROMACS coordinate file
{% endhighlight %}  

### Fix names for Amber14FF

If your system was prepared using Amber14FF, some atom-specific names are not recognized by acpype, so we need to modify the topology file:

{% highlight git %}
sed 's/ 2C / CT /g' system_solv_GMX.top | sed 's/ 3C / CT /g' | sed 's/ IP / Na+/' | sed 's/ IM / Cl-/' > system_solv_GMX_corr.top  
{% endhighlight %}

### Reorder ions.

Finally, if you added Na+ and Cl- ions in Amber, these can be in a random order in the system_solv_GMX.gro file. Check your .gro file to see if Na+ and Cl- are out of order. GROMACS needs Na+ and Cl- to all be grouped together in the coordinate file. If you need to reorder the ions, you can use the script, [reorder.py](https://github.com/jamesmccarty/amber2gromacs/blob/main/reorder.py):

{% highlight git %}
python reorder.py -i system_solv_GMX.gro -o system_solve_GMX_reorder.gro
{% endhighlight %}

Now you have a GROMACS topology and coordinate file that can be used to run MD simulations in GROMACS.
