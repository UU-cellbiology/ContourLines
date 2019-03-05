package fiji.plugin.ContourLines;

import fiji.plugin.ContourLines.Point;
import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.Toolbar;
import ij.plugin.LutLoader;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.LUT;


public class ContourLines_ implements PlugIn {
	/** main image window */
	public ImagePlus imp; 	
	/** currently active processor */
	public ImageProcessor ip;
	/** smoothing radius used to calculate gradient **/
	float dSmoothR;
	/** derivatives of gaussian (first order in x and y)**/
	ImageProcessor fx,fy;
	/**normalized vector field of gradient**/
	float [][][] fNormalizedField;
	/** integration timestep **/
	float dTimeStep;
	
	
	/** grid of occupied spaces **/
	List[][] Grid;// = new List<String>[2];
	
	/** grid size **/
	int nWG,nHG;

	/** new seed points coordinates **/
	ArrayList<Float []> fSeedPoints = new ArrayList<Float []>(); 
	
	/** separation distance betweel contour lines **/
	float dSep;
	float dGridSize;
	/** image size **/
	int nW,nH;

	/** chosen LUT name**/
	String sLUTName;
	
	int [][] RGBLutTable;
	
	/** whether to use single color or LUT **/
	boolean bSingleColor;
	
	Color systemColor;

	/** min and max intensities **/
	float fInMin;
	float fInRange;
	float fInMax;
	
	
	@Override
	public void run(String arg) {
		
		long startTime;
		long contourTime=0;
		
		int i,j;
		int nLines=0;

		Overlay image_overlay; 
		
		imp = IJ.getImage();
		if(null == imp)
		{
			IJ.noImage();
			
			return;
		}
		else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 && imp.getType() != ImagePlus.GRAY32) 
		{
		    IJ.error("8, 16 or 32-bit grayscale image  required");
		    return;
		}	
		//ask for parameters
		if(!showParametersDialog())
			return;
		
		IJ.log("Countour Lines v.0.0.2 plugin");
		//let's start measuring time		
		startTime = System.nanoTime();
				
		nW=imp.getWidth();
		nH=imp.getHeight();
		
		
		ip=imp.getProcessor();
		//ip=(FloatProcessor) ip.duplicate().convertToFloat();
		dGridSize=0.5f*dSep;
		nWG=(int)Math.ceil(nW/(float)dGridSize);
		nHG=(int)Math.ceil(nH/(float)dGridSize);

		Grid = new ArrayList[nWG][nHG];
		//fSeedPointsCoords = new float [nWG][nHG][2];
		double[][] kernel = computeKernelGaussian2D_du(dSmoothR,dSmoothR, 0.0f);		
		fx= convolve(ip, kernel);
		//new ImagePlus("dx", fx).show();
		kernel = computeKernelGaussian2D_dv(dSmoothR,dSmoothR, 0.0f);
		fy= convolve(ip, kernel);		
		//new ImagePlus("dy", fy).show();
		fNormalizedField  = new float [nW][nH][2];
		fillNormalizedField();
		image_overlay = new Overlay();
		if (bSingleColor)
			systemColor = Toolbar.getForegroundColor();
		else
		{
		    /*LUTService ls = new DefaultLUTService();
		    try {
				lut = ls.loadLUT(ls.findLUTs().get((sLUTName+".lut")));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}*/
			getRGBLutTable();
			fInMin = (float) ip.getMin();
			fInMax=(float) ip.getMax();
			fInRange = fInMax -fInMin;
			
		}
		
		
		
		image_overlay.add(buildLine((float)(nW*0.5),(float)(nH*0.5)));
		nLines++;
		//image_overlay.add(buildLine((float)(34),(float)(22)));
		imp.setOverlay(image_overlay);
		imp.updateAndRepaintWindow();
		imp.show();
		boolean goOn=true;
		float newx,newy;
		Float [] seed_point;

		while (goOn)
		{
			newx=-1;
			newy=-1;
			if(fSeedPoints.size()==0)
			{
				goOn=false;
			}
			else
			{
				boolean bFoundSeed = false;
				i=0;
				while (!bFoundSeed)
				{
					if(fSeedPoints.size()==0)
					{
						goOn = false;
						break;
					}
					seed_point= fSeedPoints.get(i);
					if(isBusy(seed_point[0],seed_point[1],0.9f*dSep))
					{
						//kick out the seed
						fSeedPoints.remove(0);
					}
					else
					{
						bFoundSeed=true;
						newx=seed_point[0];
						newy=seed_point[1];
						fSeedPoints.remove(0);
					}
				}
	
				if(goOn && bFoundSeed)
				{
					
					image_overlay.add(buildLine(newx,newy));
					nLines++;
					imp.setOverlay(image_overlay);
					imp.updateAndRepaintWindow();
					imp.show();
					IJ.log("another line "+ Integer.toString(nLines));
				}
			}
		}
		IJ.log("Done");
		contourTime = System.nanoTime() - startTime;
		IJ.log("Total time: " + String.format("%.2f",((double)Math.abs(contourTime))*0.000000001) + " s");
		
	}
	
	public PolygonRoi buildLine(float xstart, float ystart)
	{
		int [][] ownGrid;
		PolygonRoi polyline;
		float xcurr, ycurr;
		float xvel,yvel;
		ArrayList<Float> px = new ArrayList<Float>(); 
		ArrayList<Float> py = new ArrayList<Float>();
		
		ownGrid = new int [nW][nH];
		
		int gx, gy;
		int gxn, gyn;
		boolean bEnd;
		
		px.add(xstart);
		py.add(ystart);
		int nDirection;
		
		nDirection =1;
		
		for(nDirection=1;nDirection>-2;nDirection-=2)
		{
		
			//forward growth
			bEnd = false;
			xcurr=xstart;
			ycurr=ystart;
			gx = (int) Math.floor(xstart);
			gy = (int) Math.floor(ystart);
			ownGrid[gx][gy]=1;
			//addSeedPoint(xstart,ystart);
			while (!bEnd)
			{
				xvel = ((float)nDirection)*fNormalizedField[Math.round(xcurr)][Math.round(ycurr)][0];
				yvel = ((float)nDirection)*fNormalizedField[Math.round(xcurr)][Math.round(ycurr)][1];
				if(Float.isNaN((Float)xvel))
				{
					bEnd=true;
					ownGrid[gx][gy]=2;
				}
				else
				{
				if(Float.isNaN((Float)xcurr))
				{
					bEnd = true;
					ownGrid[gx][gy]=2;
				}
				else
				{
					
						
					xcurr = xcurr +xvel * dTimeStep;
					ycurr = ycurr +yvel * dTimeStep;
					//if(xcurr<0 || ycurr<0 || xcurr>(nW-1)|| ycurr>(nH-1))
					if(isBusy(xcurr,ycurr,((float)dSep)*0.5f))
					{
						bEnd = true;
						ownGrid[gx][gy]=2;
					}
					else
					{
						gxn = (int) Math.floor(xcurr);
						gyn = (int) Math.floor(ycurr);
						if(gxn<0 || gyn<0 || gxn>(nW-1)|| gyn>(nH-1))
						{
							bEnd = true;
						}
						else
						{
							if(ownGrid[gxn][gyn]==2)
							{
								bEnd = true;
							}
							else
							{
								if(gxn!=gx || gyn!=gy)
								{
									ownGrid[gx][gy]=2;
									gx=gxn;
									gy=gyn;
								}
								//addSeedPoint(xcurr,ycurr);
								if(nDirection==1)
								{
									px.add(xcurr);
									py.add(ycurr);
									
								}
								else
								{
									px.add(0,xcurr);
									py.add(0,ycurr);								
								}
		
							}
						}
					}
				}
				}
					
			}
			ownGrid[gx][gy]=2;
		}
	
		float[] floatX = new float[px.size()];
		float[] floatY = new float[px.size()];
		float nAverVal=0;
		int i = 0;
		for (i=0;i<px.size();i++) 
		{
			floatX[i]=px.get(i);
			floatY[i]=py.get(i);
			addSeedPoint(floatX[i],floatY[i]);
			markOccupied(px.get(i), py.get(i));
		}
		
		 
		
		 polyline = new PolygonRoi(floatX, floatY, Roi.POLYLINE);
		 if(bSingleColor)
		 {
			 polyline.setStrokeColor(systemColor);
		 }
		 else
		 {
			//calculating average intensity per line
			nAverVal=0;
			for (i=0;i<px.size();i++) 
			{
				nAverVal += ip.getf((int)Math.floor(floatX[i]), (int)Math.floor(floatY[i]));
			}
			nAverVal/=px.size();	
			if(nAverVal>fInMax)
				nAverVal=255;
			else
			{
				if(nAverVal<fInMin)
					nAverVal=0;
				else
						nAverVal=Math.round(255.0f*(nAverVal-fInMin)/fInRange);
			}
			
			polyline.setStrokeColor(new Color(RGBLutTable[(int)nAverVal][0],RGBLutTable[(int)nAverVal][1],RGBLutTable[(int)nAverVal][2]));
		 }
		 
		 
		
		return polyline;
		
		
	}
	void markOccupied(float xin, float yin)
	{
		int gx,gy;
		ArrayList<Point> currarr;		
		
		//old
		//bIsOccupied[Math.round(xin)][Math.round(yin)]=1;
		
		gx = (int) Math.floor(xin/dGridSize);
		gy = (int) Math.floor(yin/dGridSize);
		currarr =(ArrayList<Point>) Grid[gx][gy];
		if(currarr==null)
		{
			Grid[gx][gy] = new ArrayList<Point>();
			currarr =(ArrayList<Point>) Grid[gx][gy];
		}
		currarr.add(new Point(xin,yin));
	}
	void addSeedPoint(float xs, float ys)
	{
		float xn,yn;
		float xvel, yvel;
		//int gx, gy;
		Float[] spoint = new Float [2];
				
		xvel = fNormalizedField[Math.round(xs)][Math.round(ys)][1];
		yvel = -fNormalizedField[Math.round(xs)][Math.round(ys)][0];
		if(Float.isNaN((Float)xvel) ||Float.isNaN((Float)yvel))
		{
			return;
		}
		xn=xs+xvel*dSep*1.1f;
		yn=ys+yvel*dSep*1.1f;

		if(xn<0 || yn<0 || xn>(nW-1)|| yn>(nH-1))
		{
			//return;
		}
		else
		{
			
			spoint[0]=xn;
			spoint[1]=yn;
			fSeedPoints.add(spoint);
		}
		

		xn=xs-xvel*dSep*1.1f;
		yn=ys-yvel*dSep*1.1f;

		if(xn<0 || yn<0 || xn>(nW-1)|| yn>(nH-1))
		{
			return;
		}
		else
		{
			spoint = new Float [2];
			spoint[0]=xn;
			spoint[1]=yn;
			fSeedPoints.add(spoint);
		}
		
	}
	public boolean isBusy(float xc, float yc, float dDist)
	{
		
		int gx, gy;
		int x,y;
		int i,j,k;

		ArrayList<Point> currarr;
		Point point;
		int nRange=(int) Math.ceil(dDist/dGridSize);
		
		x= Math.round(xc);
		y= Math.round(yc);
		if(x<0 || y<0 || x>(nW-1)|| y>(nH-1))
		{
			return true;
		}
		
		gx = (int) Math.floor(xc/dGridSize);
		gy = (int) Math.floor(yc/dGridSize);
		
		for (i=gx-nRange;i<=gx+nRange;i++)
			for (j=gy-nRange;j<=gy+nRange;j++)
			{
				if(i<0 || j<0 || i>(nWG-1)|| j>(nHG-1))
				{
								
				}
				else
				{
					currarr =(ArrayList<Point>) Grid[i][j];
					if(currarr==null)
					{
					
					}
					else
					{
						for (k=0;k<currarr.size();k++)
						{
							point= currarr.get(k);
							if(point.distance(xc, yc)<dDist)
								return true;
						}
					}											
				
				}
			}
	
		return false;
		
	}
	public void fillNormalizedField()
	{
		int i,j;
		float dx,dy,len;
		
		for (i=0;i<nW;i++)
			for (j=0;j<nH;j++)
			{
				dx = fx.getf(i, j);
				dy = fy.getf(i, j);
				len = (float) Math.sqrt(dx*dx+dy*dy);
				if(len>0.00000001f)
				{
					//fNormalizedField[i][j][0]=dx/len;
					//fNormalizedField[i][j][1]=dy/len;
					fNormalizedField[i][j][0]=dy/len;
					fNormalizedField[i][j][1]=-dx/len;

				}
				else
				{
					fNormalizedField[i][j][0]=Float.NaN;
					fNormalizedField[i][j][1]=Float.NaN;					
				}
			}
		
	}
	/** Dialog with linking parameters **/
	public boolean showParametersDialog()
	{
		int nLutChoice;
		GenericDialog contourlinesD = new GenericDialog("Rendering parameters");
		String [] luts = IJ.getLuts();
		
		
		//linkingD.addChoice("Linking sequence: ",sLinkType, Prefs.get("CurveTrace.nLinkingType", "Incremental"));
		contourlinesD.addNumericField("Smoothing radius:", Prefs.get("ContourLines.dSmoothR", 2.0), 1, 3,"pixels");
		contourlinesD.addNumericField("Time Step:", Prefs.get("ContourLines.dTimeStep", 0.01), 2, 3,"s");
		contourlinesD.addNumericField("Separate distance:", Prefs.get("ContourLines.dSep", 3), 1, 3,"pixels");
		contourlinesD.addCheckbox("Use single color (current)?", Prefs.get("ContourLines.bSingleColor", true));
		contourlinesD.addChoice("Build LUT coded contours:",luts,Prefs.get("ContourLines.sLutChoice","Fire"));
		contourlinesD.setResizable(false);
		contourlinesD.showDialog();	
		if (contourlinesD.wasCanceled())
            return false;

		dSmoothR = (float) contourlinesD.getNextNumber();
		Prefs.set("ContourLines.dSmoothR", dSmoothR);
		dTimeStep = (float) contourlinesD.getNextNumber();
		Prefs.set("ContourLines.dTimeStep", dTimeStep);
		dSep = (float) contourlinesD.getNextNumber();
		Prefs.set("ContourLines.dSep", dSep);
		bSingleColor = contourlinesD.getNextBoolean();
		Prefs.set("ContourLines.bSingleColor", bSingleColor);
		nLutChoice = contourlinesD.getNextChoiceIndex();
		Prefs.set("ContourLines.sLutChoice", luts[nLutChoice]);
		sLUTName = luts[nLutChoice];

		return true;
	}
	
	public static double[][] computeKernelGaussian2D_du(double sigma_x, double sigma_y, float theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*Math.max(sigma_x, sigma_y)); // RSLV: use floor instead of round?
		int kernel_size = 1+2*kernel_radius;
		
		// compute kernel
		double[][] kernel = new double[kernel_size][kernel_size];
		for(int ky = 0; ky < kernel_size; ++ky)
		{
			int y = ky - kernel_radius;
			for(int kx = 0; kx < kernel_size; ++kx)
			{
				int x = kx - kernel_radius;
				double u = x * Math.cos(theta) - y * Math.sin(theta);
				double v = x * Math.sin(theta) + y * Math.cos(theta);
				kernel[kx][ky] = gaussian2D_dx(u, v, sigma_x, sigma_y);
			}
		}
		
		// normalize kernel
		kernel = normalize_kernel(kernel);
		
		// return kernel
		return kernel;
	}
	
	public static double[][] computeKernelGaussian2D_dv(double sigma_x, double sigma_y, float theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*Math.max(sigma_x, sigma_y)); // RSLV: use floor instead of round?
		int kernel_size = 1+2*kernel_radius;
		
		// compute kernel
		double[][] kernel = new double[kernel_size][kernel_size];
		for(int ky = 0; ky < kernel_size; ++ky)
		{
			int y = ky - kernel_radius;
			for(int kx = 0; kx < kernel_size; ++kx)
			{
				int x = kx - kernel_radius;
				double u = x * Math.cos(theta) - y * Math.sin(theta);
				double v = x * Math.sin(theta) + y * Math.cos(theta);
				kernel[kx][ky] = gaussian2D_dy(u, v, sigma_x, sigma_y);
			}
		}
		
		// normalize kernel
		kernel = normalize_kernel(kernel);
		
		// return kernel
		return kernel;
	}
	
	public static double gaussian2D_dx(double x, double y, double sigma_x, double sigma_y)
	{
		return ((-x)/(2*Math.PI*Math.pow(sigma_x, 3)*sigma_y))*Math.exp(-0.5*((x*x)/(sigma_x*sigma_x)+(y*y)/(sigma_y*sigma_y)));
	}
	
	public static double gaussian2D_dy(double x, double y, double sigma_x, double sigma_y)
	{
		//return gaussian2D_dx(y, x, sigma); // NOTE x and y are swapped
		return ((-y)/(2*Math.PI*sigma_x*Math.pow(sigma_y, 3)))*Math.exp(-0.5*((x*x)/(sigma_x*sigma_x)+(y*y)/(sigma_y*sigma_y)));
	}
	public static double[][] normalize_kernel(double[][] kernel)
	{
		// calculate sum of components
		double sum = 0.0;
		for(int kx = 0; kx < kernel.length; ++kx)
		{
			for(int ky = 0; ky < kernel[kx].length; ++ky)
			{
				sum += Math.abs(kernel[kx][ky]); // NOTE: use abs to normalize symmetrical kernel with a positive and negative lobe
			}
		}
		
		// avoid division by zero
		if(sum == 0.0) { return kernel; }
		
		// calculate scale factor
		double scale_factor = 1 / sum;
		
		// scale components
		for(int kx = 0; kx < kernel.length; ++kx)
		{
			for(int ky = 0; ky < kernel[kx].length; ++ky)
			{
				kernel[kx][ky] *= scale_factor;
			}
		}
		
		// return normalized kernel
//		ImageProcessor kernel_ip = kernelToImage(kernel);
//		ImagePlus kernel_imp = new ImagePlus("normalized kernel", kernel_ip);
//		kernel_imp.resetDisplayRange();
//		kernel_imp.show();
		return kernel;
	}
	
	public static ImageProcessor convolve(ImageProcessor ip, double[][] kernel)
	{
		// get image and kernel sizes
		int image_width = ip.getWidth();
		int image_height = ip.getHeight();
		
		int kernel_width = kernel.length;
		int kernel_height = kernel_width; // NOTE: assume square kernel
		int kernel_half_width = (int)Math.floor(0.5 * kernel_width);
		int kernel_half_height = kernel_half_width;
		
		// convert input image processor to float
		ImageProcessor ip_inp = ip.convertToFloat(); // RSLV: duplicate first?
		
		// create new empty output float processor
		ImageProcessor ip_res = new FloatProcessor(image_width, image_height);
		
		// convolve input image with kernel
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				double kernel_product = 0.0;
				for(int ky = 0; ky < kernel_height; ++ky)
				{
					int ppy = py + ky - kernel_half_height;
					if(ppy < 0) ppy = 0; // clamp at border
					if(ppy >= image_height) ppy = image_height - 1; // clamp at border
					for(int kx = 0; kx < kernel_width; ++kx)
					{
						int ppx = px + kx - kernel_half_width;
						if(ppx < 0) ppx = 0; // clamp at border
						if(ppx >= image_width) ppx = image_width - 1; // clamp at border
						kernel_product += ip_inp.getf(ppx, ppy) * kernel[kx][ky];
					}
				}
				ip_res.setf(px, py, (float)kernel_product);
			}
		}
		
		// return 
		ip_res.resetMinAndMax();
		return ip_res;
	}
	
	/** function gets LUT specified by sZLUTName in settings
	 * and returns 256x3 table map in HSB format */
	void  getRGBLutTable()
	{
		int i,j;
		
		
		
		int [] onepix; 
		RGBLutTable = new int[256][3];
		ByteProcessor ish = new ByteProcessor(256,10);
		for (i=0; i<256; i++)
			for (j=0; j<10; j++)
				ish.putPixel(i, j, i);
		ImagePlus ccc = new ImagePlus("test",ish);
		ccc.show();
		IJ.run(sLUTName);
		IJ.run("RGB Color");
		//ipLUT= (ColorProcessor) ccc.getProcessor();
		ccc.setSlice(1);
		for(i=0;i<256;i++)
		{
			
			onepix= ccc.getPixel(i, 2);
			//rgbtable[i]=ccc.getPixel(i, 1);
			RGBLutTable[i][0]=onepix[0];
			RGBLutTable[i][1]=onepix[1];
			RGBLutTable[i][2]=onepix[2];
		}
		ccc.changes=false;
		ccc.close();
		
		return;
	}

}
