# ContourLines
[ImageJ](https://imagej.nih.gov/ij/)/[FIJI](http://fiji.sc/) plugin generating contour lines with equal spacing on top of an image (using overlay).  
Based/inspired by [streamlines](https://github.com/anvaka/streamlines) project by [anvaka](https://github.com/anvaka) and [this publication](http://web.cs.ucdavis.edu/~ma/SIGGRAPH02/course23/notes/papers/Jobard.pdf) (use it for citation).
   
![contourexample](http://katpyxa.info/software/ContourLines/CL2.gif "logo")

## How to install plugin

You need to download and install [ImageJ](https://imagej.nih.gov/ij/download.html) or [FIJI](http://fiji.sc/#download) on your computer first.

For FIJI
* add https://sites.imagej.net/Ekatrukha/ to the list of update sites, as follows
* go to Help -> Update and press "Manage update sites" button
* press "Add update site" button and put the following link there https://sites.imagej.net/Ekatrukha/

For ImageJ
* Download and copy [the latest version of plugin](https://github.com/ekatrukha/ContourLines/blob/master/target/ContourLines_-0.0.4.jar?raw=true) into the *plugins* folder of your ImageJ. (for example, in Windows look for *C:\Program Files\ImageJ\plugins*)

Restart ImageJ/Fiji and plugin will appear in *Plugins* menu

## Parameters

After plugin launch, rendering parameters window will appear  

![paramwindows](http://katpyxa.info/software/ContourLines/CL_parameters_dialog_v.0.0.4.png "parameters window")

### Smoothing radius
To find the "vector field of an image", plugin calculates convolution of image with [derivative of Gaussian](http://campar.in.tum.de/Chair/HaukeHeibelGaussianDerivatives) in *x* and *y* to estimate intensity gradient at each point. You can specify smoothing radius (SD of Gaussian). The larger the value, less small intensity details are available. For example, here are two images, the left one is build with smoothing radius of 2.0, while the right one with smoothing radius of 10.0 pixels.

![Rad2](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_x.png "Rad2")  ![Rad10](http://katpyxa.info/software/ContourLines/smoothing_10_line_0.05_distance_3_x.png "Rad10"). 

Notice that the small spot in the right top corner "disappeared".

### Line integration step

This parameter defines, how precisely contour generation routine follows gradient vector field of the image. Smaller value result in more precise contour lines, but also take longer time to generate them. Example below illustrates integration step of 0.05 (left) versus 0.5 pixels (right).

![Step005](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_x.png "Step005")  ![Step05](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.5_distance_3_x.png "Step05"). 

### Distance between lines

Well, how densely you want contour lines to be plotted. Left image is 3 px, right image is 10 px.

![dist3](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_x.png "Step005")  ![dist10](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_10_x.png "dist10"). 

### Stop line at fraction of distance

This parameter is a fraction of previous and it defines when the line integration will stop. I.e. end of each line cannot come closer than this distance to already existing lines. For example, on the left image parameter's value is 0.5, while on the right it is 0.1.  

![stop05](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_x.png "Stop05")  ![stop01](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_single_color_end_0.1.png "stop01"). 

### Use single color

With this option checked, plugin will use current selected color (in toolbar or Color Picker) to build lines. Example:

![singlecolor](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_single_color_x.png "single color")

### Use color lookup tables (LUT)

Plugin will use one of the installed ImageJ/Fiji lookup tables (LUTs) to color code each contour with its average intensity. Pay attention, that current minimum and maximum of intensity would be taken from "Brightness/Contrast" settings of the image. In works better if your image uses "dark-to-bright" LUT and you check "Invert LUT" option for contours. In this case brighter lines would be on top of dark image areas and vice versa. In the example below the left image uses inverted "Thermal" LUT for contours, while contour LUT of the image on the right is not inverted.

![invertlut](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_x.png "invertlut")  ![notinvertlut](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_not_inverted_thermal_x.png "notinvertlut"). 

### Remove open contours

Plugin will remove open contours (contours with free ends). Examples are below, left image includes open contours, while the right one not.

![withopen](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_x.png "withopen")  ![noopen](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_single_color_closed_only.png "noopen"). 

## Updates history

* v.0.0.4 added criteria for line integration end. Added a possibility to remove open contours.

---
Developed in [Cell Biology group](http://cellbiology.science.uu.nl/) of Utrecht University.  
Email katpyxa @ gmail.com for any questions/comments/suggestions.
