mesh_segmentation
=================

A simple python addon for blender using spectral clustering to segment meshes. Uses [numpy](http://www.numpy.org/) and [scipy](http://www.scipy.org/) for matrix calculations.
Developed as a project for the spectral-clustering course of University Bremen. Based on the paper "Segmentation of 3D Meshes through Spectral Clustering" by Rong Liu and Hao Zhang.
### Usage

 - Install [scipy](http://scipy.org/install.html) in the Python installation that is used by Blender. On Linux this should usually be the system Python, so simply installing *scipy* with your package manager should work. On Windows there will usually be a separate Python installation coming with Blender, so it's a bit more tricky. Perhaps the simplest way is to open a command prompt or power shell with admin rights, navigate to "YOUR_BLENDER_PATH/python/bin" and run `./python -m ensurepip` followed by `./python -m pip install scipy`
 - Put the "mesh_segmentation" folder inside of "src" to "YOUR_BLENDER_PATH/python/scripts/addons"
 - Run Blender, go to "User Preferences -> Addons" and search for "mesh segmentation"
 - Activate the plugin by ticking the box at the right.

Now you can select an object in the 3D viewport, hit spacebar and search for "Segment Mesh". Clicking on that entry shows the addon in a popup.

<img src="example-results/hand.png" alt="Hand">
You can find more examples in the example-results folder.
