# Installing

To install the plugin, simply place the Subchondral_Thickness.java file into the plugins folder of ImageJ. 
> If you place it inside a folder within 'plugins' it will be within that submenu when running ImageJ. 

Start ImageJ, click
on the 'Plugins' menu, and select 'Compile and Run...'. 

In the resulting dialog box, navagate to where you saved Subchondral_Thickness.java
and select it. 

The plugin should run and the next time you start ImageJ, it will be an option directly from the Plugins menu.

## Using the Plugin

To use the plugin, follow these steps:

1. Import a sequence of DICOM images, via File -> Import -> Image Sequence...

2. Run the Subchondral Thickness plugin

3. On the image you want to start with, click 4 points to create a Region of Interest

	> **Note** When selecting your points, the screen may not show the most recent point you clicked, a fix is in the works, but for now be confident that it did register your click
	> 
	> **Note** When selecting your first ROI, always start with the higher image in the stack (lower number). Curently the plugin can only process from the first image down to the second. Again, a fix is in the works

4. After selecting the points, click 'Ok' in the dialog

5. On the image you want to end with, click 4 points to create a Region of Interest 

	> **NOTE** Currently, you must select the points in the same order as the first image, i.e. if you clicked the top left corner first, click the top left corner of the second region first as well. The order matters, as the plugin linearly interpolates the points
			 between matching verticies 
	
	> **Note** When selecting your points, the screen may not show the most recent point you clicked, a fix is in the works, but for now be confident that it did register your click

6. After selecting the points, click 'Ok' in the dialog

7. Input the file (including the file path) where you want the results to be saved to. If the file exists, the new results will be appended
   to the end. If the file doesn't exist, it will be created.

8. Adjust the threshold slider to your satisfaction. Checking the Preview box will show you what the threshold will look like

9. Click 'Ok'

10. On the dialog that pops up, asking if it should run on all images in the stack, click 'No'. The plugin automatically moves through the stack. (Working on surpressing this dialog)

10. The plugin will now run and display the results and save them to the file specified


