# ContourLines
[ImageJ](https://imagej.nih.gov/ij/)/[FIJI](http://fiji.sc/) plugin generating contour lines with equal spacing on top of an image (using overlay).  
Based/inspired by [streamlines](https://github.com/anvaka/streamlines) project by [anvaka](https://github.com/anvaka) and [this publication](http://web.cs.ucdavis.edu/~ma/SIGGRAPH02/course23/notes/papers/Jobard.pdf) (use it for citation).
   
![contourexample](http://katpyxa.info/software/ContourLines/CL2.gif "logo")

## How to install plugin

* You need to download and install [ImageJ](https://imagej.nih.gov/ij/download.html) or [FIJI](http://fiji.sc/#download) on your computer first.
* Download and copy [the latest version of plugin](https://github.com/ekatrukha/ContourLines/blob/master/target/ContourLines_-0.0.3.jar?raw=true) into the plugins folder of your ImageJ of Fiji. (for example, in Windows look for *C:\Program Files\ImageJ\plugins*)
* restart ImageJ/Fiji
* plugin will appear in *Plugins* menu

## Parameters

After plugin launch, rendering parameters window will appear  

![paramwindows](http://katpyxa.info/software/ContourLines/CL_parameters_dialog.png "parameters window")

### Smoothing radius
To find the "vector field of an image", plugin calculates convolution of image with [derivative of Gaussian](http://campar.in.tum.de/Chair/HaukeHeibelGaussianDerivatives) in *x* and *y* to estimate intensity gradient at each point. You can specify smoothing radius (SD of Gaussian). The larger the value, less small intensity details are available. For example, here are two images, the left one is build with smoothing radius of 2.0, while the right one with smoothing radius of 10.0 pixels.

![Rad2](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3.png "Rad2")  ![Rad10](http://katpyxa.info/software/ContourLines/smoothing_10_line_0.05_distance_3.png "Rad10"). 

Notice that the small spot in the right top corner "disappeared".

### Line integration step

This parameter defines, how precisely contour generation routine follows gradient vector field of the image. Smaller value result in more precise contour lines, but also take longer time to generate them. Example below illustrates integration step of 0.05 (left) versus 0.5 pixels (right).

![Step005](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3.png "Step005")  ![Step05](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.5_distance_3.png "Step05"). 

### Distance between lines

Well, how densely you want contour lines to be plotted. Left image is 3 px, right image is 10 px.

![dist3](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3.png "Step005")  ![dist10](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_10.png "dist10"). 

### Use single color

With this option checked, plugin will use current selected color (in toolbar or Color Picker) to build lines. Example:

![singlecolor](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_single_color.png "single color")

### Use color lookup tables (LUT)

Plugin will use one of the installed ImageJ/Fiji lookup tables (LUTs) to color code each contour with its average intensity. Pay attention, that current minimum and maximum of intensity would be taken from "Brightness/Contrast" settings of the image. In works better if your image uses "dark-to-bright" LUT and you check "Invert LUT" option for contours. In this case brighter lines would be on top of dark image areas and vice versa. In the example below the left image uses inverted "Thermal" LUT for contours, while contour LUT of the image on the right is not inverted.

![invertlut](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3.png "invertlut")  ![notinvertlut](http://katpyxa.info/software/ContourLines/smoothing_2_line_0.05_distance_3_not_inverted_thermal.png "notinvertlut"). 

## License
MIT

---
Developed in [Cell Biology group](http://cellbiology.science.uu.nl/) of Utrecht University.  
Email katpyxa @ gmail.com for any questions/comments/suggestions.