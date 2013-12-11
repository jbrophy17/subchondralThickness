import ij.*;
import ij.gui.*;
import ij.measure.Calibration;
import ij.process.*;
import ij.plugin.filter.*;

import java.awt.*;
import java.util.*;
import java.awt.event.*;
import java.io.FileWriter;
import java.io.IOException;

/**
 *	Plugin to measure the average thickness of bone over a specified
 *  3D region of interest. Measures the average 2D thickness on each
 *  slice of a stack of DICOM images and then averages all the slices.
 *  Upon completion, displays result and outputs it to the specified
 *  output file.
 *
 *	@author	John Brophy
 *	@author	Washington University in St. Louis
 */


public class Subchondral_Thickness implements ExtendedPlugInFilter, DialogListener, MouseListener {
	ImagePlus imp = null;
	int firstROIindex = -1;
	Roi firstROI = null;
	int secondROIindex = -1;
	Roi secondROI = null;
	private static int intThresh = 127;
	private static int intThreshMin = 1;
	private static int intThreshMax = 1;
	private static int minRange = 0;
	private static int maxRange = 1;
	
	private static String filepath = "";
	private static String title = "";
	
	private int invalidSlices = 0;
	
	boolean isPreview;
	GenericDialog thresholdDialog;
	
	Overlay roiOverlay;
	
	double physicalX; 
	double physicalY; 
	double physicalZ; 
	String spaceUnit; 
	
	ImageCanvas canvas;
	
	//Booleans required because accessing the already created selection tools is basicaly impossible
	boolean isFirstRoi = false;
	boolean isSecondRoi = false;
	float[] firstXPoints = new float[4];
	float[] firstYPoints = new float[4];
	float[] secondXPoints = new float[4];
	float[] secondYPoints = new float[4];
	int counter = 0;

	public int setup(String arg, ImagePlus imp) {
		if (imp == null) {
			IJ.noImage();
			return DONE;
		}
		if (arg.equals("about")) {
			showAbout();
			return DONE;
		}
		//Get the stack of images
		ImageStack curStack = imp.getImageStack();
		//Make sure its a valid stack
		if (curStack.getSize() < 2) {
			showImageStackError();
			return DONE;
		}
		
		title = imp.getTitle();
		IJ.register(Subchondral_Thickness.class);
		ImageWindow win = imp.getWindow();
		canvas = win.getCanvas();
		canvas.addMouseListener(this);
		
		roiOverlay = new Overlay();
		imp.setOverlay(roiOverlay);
		while (firstROIindex == -1 || firstROI == null) {
			isFirstRoi = true;
			new WaitForUserDialog("Select the first ROI").show();
			if(isFirstRoi) {
				IJ.showMessage("You need to specify 4 points");
				counter = 0;
				roiOverlay.clear();
				continue;
			}
			roiOverlay.clear();
			firstROIindex = imp.getCurrentSlice();
				
			firstROI = new PolygonRoi(firstXPoints,firstYPoints,4,Roi.FREEROI);
			
		}

		while (secondROIindex == -1 || secondROI == null || secondROI == firstROI) {
			isSecondRoi = true;
			new WaitForUserDialog("Select the second ROI").show();
			if(isSecondRoi) {
				IJ.showMessage("You need to specify 4 points");
				counter = 0;
				roiOverlay.clear();
				continue;
			}
			roiOverlay.clear();
			secondROIindex = imp.getCurrentSlice();
				
			secondROI = new PolygonRoi(secondXPoints,secondYPoints,4,Roi.FREEROI);
			
		}
		
		this.imp = imp;
		Calibration cal = imp.getCalibration();
		double[] coeficients = cal.getCoefficients();
		intThreshMin = (int) imp.getDisplayRangeMin();// + 500;
		intThreshMax = (int) imp.getDisplayRangeMax();// - 100;
		maxRange = intThreshMax;
		minRange = intThreshMin;
		
		intThreshMin += coeficients[0];
		intThreshMax += coeficients[0];
		
		intThreshMin = (int) (intThreshMin / coeficients[1]);
		intThreshMax = (int) (intThreshMax / coeficients[1]);
		
		intThresh = (intThreshMin + intThreshMax) / 2;
		
		
		physicalX = cal.pixelWidth; 
		physicalY = cal.pixelHeight; 
		physicalZ = cal.pixelDepth; 
		spaceUnit = cal.getUnit(); 
		
		return DOES_16 + STACK_REQUIRED + SUPPORTS_MASKING;
	}
	
	public void run(ImageProcessor arg) {
		if (thresholdDialog.wasOKed()) {
			isPreview = false;
		}
		invalidSlices = 0;
		
		ImageStack curStack = imp.getImageStack();
		float[] xOffset = new float[4];
		float[] yOffset = new float[4];
		int size = secondROIindex - firstROIindex + 1;
		
		xOffset[0] = (firstXPoints[0] - secondXPoints[0] ) / size;
		yOffset[0] = (firstYPoints[0] - secondYPoints[0] ) / size;
		
		xOffset[1] = (firstXPoints[1] - secondXPoints[1] ) / size;
		yOffset[1] = (firstYPoints[1] - secondYPoints[1] ) / size;
		
		xOffset[2] = (firstXPoints[2] - secondXPoints[2] ) / size;
		yOffset[2] = (firstYPoints[2] - secondYPoints[2] ) / size;
		
		xOffset[3] = (firstXPoints[3] - secondXPoints[3] ) / size;
		yOffset[3] = (firstYPoints[3] - secondYPoints[3] ) / size;
		
		double totalThickness = 0.0;
		
		for( int i = 0; i < size; ++i) {
			
			if(isPreview) {
				if (i != size - 1) {
					continue;
				}
			}
			
			ImageProcessor holdIp = curStack.getProcessor(i + firstROIindex);
			float[] holdXCoords = new float[4];
			float[] holdYCoords = new float[4];
			
			for ( int j = 0; j < 4; ++j) {
				holdXCoords[j] = firstXPoints[j] - i*xOffset[j];
				holdYCoords[j] = firstYPoints[j] - i*yOffset[j];
			}
			PolygonRoi holdRoi = new PolygonRoi(holdXCoords,holdYCoords,4,Roi.FREEROI);
			byte[] holdBinaryArray = thresholdImage(holdIp, holdRoi);
			
			ByteProcessor holdByteProcessor = new ByteProcessor(holdRoi.getBounds().width, holdRoi.getBounds().height, holdBinaryArray);
			BinaryProcessor holdSkeletonBinaryProcessor = new BinaryProcessor((ByteProcessor)holdByteProcessor.duplicate());
			BinaryProcessor holdOutlineBinaryProcessor = new BinaryProcessor((ByteProcessor)holdByteProcessor.duplicate());

			holdSkeletonBinaryProcessor.invert();
			holdSkeletonBinaryProcessor.skeletonize();
			holdOutlineBinaryProcessor.invert();
			holdOutlineBinaryProcessor.outline();

			double holdDistance = findDistanceBetweenEdges(holdOutlineBinaryProcessor, holdRoi);
			if (holdDistance < 0) {
				if (!isPreview) {
					IJ.log("An error has occured while processing the image at slice " + i + firstROIindex + ". \n"
						+ "This could be due to your threshold being too small / too big or \n" + 
						"your ROI not having enough padding between its edge and the shell of the bone. \n" +
						"Please Note that this slice has been excluded from the analysis");
				}
				continue;
			}
			totalThickness += holdDistance;
		}
		
		size -= invalidSlices;
		
		totalThickness = totalThickness / size;
		
		double averageThickness = totalThickness * physicalX;
		double depth = size * physicalZ;
		
		if (!isPreview) {
			FileWriter fw = null;
			try
			{
			    if ( filepath.equals("")) {
			    	IJ.log("No filename specified, unable to save results");
			    } else {
			    	fw = new FileWriter(filepath,true); //the true will append the new data to the end of the existing file
			    	fw.write("Image: " + title + " " + "Average thickness: " + averageThickness+ " " + spaceUnit + " depth: " + depth + " " + spaceUnit+"\n");
			    	
			    }
			}
			catch(IOException ioe)
			{
			    System.err.println("IOException: " + ioe.getMessage());
			} finally {
				if( fw != null) {
					try {
						fw.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
				
			}
			
			IJ.showMessage("The Average Thickness of your sample is " + averageThickness + " " + spaceUnit + " \n" +
				"The depth of your sample is " + depth + " " + spaceUnit);
		}
		
		
	}
	
	
	byte[] thresholdImage(ImageProcessor ip, Roi roi) {
		short[] pixels = (short[])ip.getPixels();
		
		if (roi!=null && !roi.isArea()) roi = null;
		
		ImageProcessor mask = roi!=null?roi.getMask():null;
		Rectangle r = roi!=null?roi.getBounds():new Rectangle(0,0,ip.getWidth(),ip.getHeight());
		int i = 0;
		int offset = ip.getWidth();
		byte[] binaryArray = new byte[r.width * r.height];
		
		if (mask == null) {
			IJ.log("null mask");
		} 
		
		for (int y=0; y<r.height; y++) {
			i = offset*(y+r.y);
			for (int x=0; x<r.width; x++) {
				if (mask==null||mask.getPixel(x,y)!=0) {
					if (ip.getPixelValue(x+r.x, y+r.y) < intThresh) {
						//make it black
						pixels[i + x + r.x] = (short)minRange;
						binaryArray[x + y*r.width] = (byte) 0;
					} else {
						//make it white
						pixels[i + x + r.x] = (short)maxRange;
						binaryArray[x + y*r.width] = (byte) 255;
					}
				}
			}
		}
		
		return binaryArray;
	}
	
	double findDistanceBetweenEdges(BinaryProcessor outline, Roi roi) {
		Map<Integer, Set<Point>> outlineMap = processPixels(outline, roi);
		//Threshold or ROI is invalid
		if (outlineMap.keySet().size() < 2) {
			if (!isPreview) {
				invalidSlices++;
			} else {
				if (intThresh - intThreshMin < intThreshMax - intThresh  && intThresh > intThreshMin) {
					intThreshMin = intThresh + 1;
					intThresh += 1;
				} else if (intThresh - intThreshMin >= intThreshMax - intThresh  &&  intThresh < intThreshMax) {
					intThresh -= 1;
					intThreshMax = intThresh;	
				} else {
					return - 1.0;
				}
				
				//TODO: Fix this warning
				@SuppressWarnings("unchecked")
				Vector<Scrollbar> sliders = (Vector<Scrollbar>)thresholdDialog.getSliders();
				for(Scrollbar s : sliders) {
					s.setValue(intThresh);
					s.setMaximum(intThreshMax);
					s.setMinimum(intThreshMin);
				}
				thresholdDialog.getPreviewCheckbox().setState(true);
			}
			return - 1.0;
		}
		
		//Only compares the two biggest
		int maxSize = 0;
		Integer maxKey = 1;
		int secondSize = 0;
		Integer secondKey = 2;
		if (outlineMap.keySet().size() > 2) {
			for (Integer i : outlineMap.keySet()) {
				if (outlineMap.get(i).size() > maxSize) {
					secondSize = maxSize;
					secondKey = maxKey;
					maxSize = outlineMap.get(i).size();
					maxKey = i;
				} else if (outlineMap.get(i).size() > secondSize) {
					secondSize = outlineMap.get(i).size();
					secondKey = i;
				}
			}
		}
		
		double distance = 0.0;
		
		for (Point oPoint : outlineMap.get(maxKey)) {
			double minDist = 9999;
			for (Point iPoint : outlineMap.get(secondKey)) {
				double holdDist = oPoint.distanceSq(iPoint);
				
				if (holdDist < minDist) {
					minDist = holdDist;
				}
			}
			distance += Math.sqrt(minDist) + 1;
		}
		
		return distance / outlineMap.get(maxKey).size();
	}
	
	Map<Integer,Set<Point>> processPixels(BinaryProcessor bp, Roi roi) {

		byte[] bpPixels = (byte[])bp.getPixels();
		
		ImageProcessor mask = roi!=null?roi.getMask():null;
		ImageProcessor maskEdge = mask.duplicate();
		maskEdge.findEdges();
		Rectangle r = roi!=null?roi.getBounds():new Rectangle(0,0,bp.getWidth(),bp.getHeight());
		
		int offset = bp.getWidth();
		int i =0;
		
		for (int y=0; y<r.height; y++) {
			i = offset*(y);
			for (int x=0; x<r.width; x++) {
				if ( (x!=0 && y!=0 && x!=r.width-1 && y!=r.height-1) && (mask==null||(mask.getPixel(x,y)!=0 && maskEdge.getPixel(x, y)==0)) && bp.getPixelValue(x, y) == 0 ) {
					bpPixels[i+x] = 0;
					
				} else {
					bpPixels[i+x] = (byte) 255;
				}
			}
		}
		
		//Connected Components
		int currentLabel = 1;
		Map<Integer, Set<Point>> resultMap = new HashMap<Integer, Set<Point>>();
		byte[][] labels = new byte[bp.getWidth()][bp.getHeight()];
		Queue<Point> Q = new LinkedList<Point>();
		
		for (int x = 0; x < bp.getWidth(); ++x) {
			for (int y = 0; y < bp.getHeight(); ++y) {
				if (labels[x][y] == 0 && bp.getPixelValue(x, y) == 0) {
					Set<Point> resultSet = new HashSet<Point>();
					resultSet.add(new Point(x,y));
					labels[x][y] = (byte) currentLabel;
					Q.add(new Point(x,y));
					
					while (!Q.isEmpty()) {
						Point hold = Q.poll();
						for(Point p : findNeighbors(hold.x,hold.y,bp.getWidth(), bp.getHeight(),bp,labels)) {
							resultSet.add(p);
							labels[p.x][p.y]= (byte) currentLabel; 
							Q.add(p);	
						}
					}
					resultMap.put(currentLabel, resultSet);
					currentLabel++;
				}
			}	
		}
		
		return resultMap;
	}
	
	ArrayList<Point> findNeighbors(int x, int y, int width, int height, ImageProcessor ip, byte[][]labels) {
		ArrayList<Point> returnPoints = new ArrayList<Point>();
		
		for (int yi = -1; yi <= 1; ++yi) {
			for (int xi = -1; xi <= 1; ++xi) {
				if ( !(xi == 0 && yi == 0) ) {
					//makes it 4 connectivity, can change easily
					if ( xi == 0 || yi == 0) {
						int xIndex = xi + x;
						int yIndex = yi + y;
						if (xIndex < 0 || yIndex < 0 || xIndex >= width || yIndex>= height) {
							continue;
						} else if (labels[xIndex][yIndex] == 0) {
							if (ip.getPixelValue(xIndex, yIndex) == 0) {
								returnPoints.add(new Point(xIndex, yIndex));
							}
						}
					}
				}
			}
		}
		
		return returnPoints;
	}

	void showAbout() {
		IJ.showMessage("About Subchondral_Thickness",
		"This Plugin is used to measure the thickness of \n" +
		"subchondral bones. To use, select two regions of \n" +
		"interest and the program will take care of the rest."
		);
	}
	
	void showImageStackError() {
		IJ.showMessage("Invalid Image Stack",
		"This Plugin requires you to specify a stack of images. \n" +
		"Please open a stack of images and try again."
		);
	}

	/*
	 * ********************
	 * Dialog Box Methods
	 * ********************
	*/
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		intThresh = (int)gd.getNextNumber();
		filepath = gd.getNextString();
		isPreview = gd.getPreviewCheckbox().getState();
		
    	return true;
    }

	public void setNPasses(int arg0) {}

	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
		GenericDialog gd = new GenericDialog(command);
		thresholdDialog = gd;
		isPreview = false;
		gd.addStringField("Location of Output File", "./output.txt");
		gd.addSlider("Threshold: ", intThreshMin, intThreshMax, intThresh);
		
		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled()) { 
			return DONE;
		}
		if (gd.wasOKed()) {
			isPreview = false;
		}
		IJ.register(this.getClass());
		return IJ.setupDialog(imp, DOES_16);
    }

	/*
	 * ********************
	 * Mouse Listener Methods
	 * ********************
	*/
	//TODO: Fix so the screen always updates on click
	@Override
	public void mouseClicked(MouseEvent e) {
		double magnification = canvas.getMagnification();
		
		if (isFirstRoi) {
			firstXPoints[counter] = (float) (e.getX() / magnification);
			firstYPoints[counter] = (float) (e.getY() / magnification);
			
			if(counter > 0) {
				Roi roi = new Line(firstXPoints[counter-1], firstYPoints[counter-1], firstXPoints[counter], firstYPoints[counter]);
				roiOverlay.add(roi);
			}
			Roi point = new OvalRoi(firstXPoints[counter] - 1, firstYPoints[counter] - 1, 3,3);
			point.setFillColor(Color.YELLOW);
			roiOverlay.add(point);
			counter++;
			if(counter > 3) {
				isFirstRoi = false;
				Roi roi = new Line(firstXPoints[counter-1], firstYPoints[counter-1], firstXPoints[0], firstYPoints[0]);
				roiOverlay.add(roi);
				counter = 0;
			}
			imp.repaintWindow();
		} else if (isSecondRoi) {
			secondXPoints[counter] = (float) (e.getX() / magnification);
			secondYPoints[counter] = (float) (e.getY() / magnification);
			if(counter > 0) {
				Roi roi = new Line(secondXPoints[counter-1], secondYPoints[counter-1], secondXPoints[counter], secondYPoints[counter]);
				roiOverlay.add(roi);
			}
			Roi point = new OvalRoi(secondXPoints[counter] - 1, secondYPoints[counter] - 1, 3,3);
			point.setFillColor(Color.YELLOW);
			roiOverlay.add(point);
			counter++;
			if(counter > 3) {
				isSecondRoi = false;
				Roi roi = new Line(secondXPoints[counter-1], secondYPoints[counter-1], secondXPoints[0], secondYPoints[0]);
				roiOverlay.add(roi);
				counter = 0;
			}
			imp.repaintWindow();
			
		}
		
	}

	@Override
	public void mouseEntered(MouseEvent arg0) {}

	@Override
	public void mouseExited(MouseEvent arg0) {}

	@Override
	public void mousePressed(MouseEvent arg0) {}

	@Override
	public void mouseReleased(MouseEvent arg0) {}

}


