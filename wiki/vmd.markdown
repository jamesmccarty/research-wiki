---
title:  "High Quality Figures in VMD"
date:   2021-11-18
layout: "git-wiki-post"
---

## Making High Quality Figures in VMD ##

### To make an image with AO (ambient occlusion lighting) and depth cueing ###

Try changing the Display Rendermode to GLSL (requires modern graphics card with OpenGL capability):

{% highlight git %}
Display --> Rendermode --> GLSL
{% endhighlight %}

Change the Display Settings:

{% highlight git %}
Display --> Display Settings ...
{% endhighlight %}

Choose:
* Shadows On
* Amb Occl On
* Cue Mode: Linear
* Cue Start: 1.75
* Cue End: 3.0

Now change the background to white and remove the axis:

{% highlight git %}
Graphics --> Colors: Display : Background : White  
{% endhighlight %}

{% highlight git %}
Display --> Axes --> Off  
{% endhighlight %}

### Experiment with different renderings ###

Try AOChalky (for tachyon ray tracer rendering). Glossy material looks nice in snapshot rendering.

### Experiment with depth cuing and orthographic perspective ###

Under the *Display* tab uncheck *Depth Cueing* and change from *Perspective* to *Orthographic*

### Render high resolution image ###

Use Tachyon module in File Render Controls:

{% highlight git %}
File --> Render ...  --> Tachyon  
{% endhighlight %}

In the same File Render Controls panel add the flag *-res xxx yyy* defining the number of pixels of the image. For example:

{% highlight git %}
"/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 1600 1200 -o %s.tga  
{% endhighlight %}
