/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import implementation.InterpolationMethods;
import java.awt.image.BufferedImage;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import util.ControlOptions;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     

    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }


    void slicer(double[] viewMatrix) { 

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];           //*view vec has the values of the voxels for a ray (pen-> bottle 1,1,0 1,1,1 1,1,2),
        double[] uVec = new double[3];              //*dimensions of plane
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2]; //* voolume center is half of every dimensions(x,y,z,)

                int val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;     //* max 3000 for example value =3000 so this is the most intense voxel
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;     //multiply with 255 for the color
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor); //* you basically "store" the values in the finalimage that you are seeing
            }
        }

    }
    
    //-------------------------------------------------------------------------------------------------------------------
    void mip(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];           //*view vec has the values of the voxels for a ray (pen-> bottle 1,1,0 1,1,1 1,1,2),
        double[] uVec = new double[3];              //*dimensions of plane
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

//        short nSlices = 50; // Default to 50 slices
        
        // My interval is upper bounded by the greatest dimension since we need
        // to respect the interval for all dimensions (assumption we make).
        List<Integer> allDimensions = new LinkedList<Integer>();
        allDimensions.add(volume.getDimX());
        allDimensions.add(volume.getDimY());
        allDimensions.add(volume.getDimZ());
        Collections.sort(allDimensions);
        
        double interval = Math.sqrt( Math.pow(allDimensions.get(2), 2) + Math.pow(allDimensions.get(1), 2) ) / ControlOptions.N_SLICES;
        
        
        int maxVox = 0;
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                maxVox = 0;
                
                for (double t = interval*(ControlOptions.N_SLICES/-2); t <= interval*(ControlOptions.N_SLICES/2) ; t=t+interval ){
                    
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + viewVec[0] * t + volumeCenter[0];
                     pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + viewVec[1] * t + volumeCenter[1] ;
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + viewVec[2] * t + + volumeCenter[2] ;
                    
                    try {
                        int val = (int)InterpolationMethods.getVoxel(pixelCoord, volume);

                        if (val > maxVox){
                            maxVox = val;
                        }
                    } catch (ArrayIndexOutOfBoundsException ex ){
                        /* If a voxel is out of bounds (meaning that we are going outside the limits, then we 
                        just ignore. This is possible because we take the maximum dimension as reference. */
                        
                    }                    
                }
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxVox/max;     //* max 3000 for example value =3000 so this is the most intense voxel
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxVox > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;     //multiply with 255 for the color
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor); //* you basically "store" the values in the finalimage that you are seeing
            }
        }

    }
    //------------------------------------------------------------------------------------------------------------------


    //-------------------------------------------------------------------------------------------------------------------
    void compositing(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];           //*view vec has the values of the voxels for a ray (pen-> bottle 1,1,0 1,1,1 1,1,2),
        double[] uVec = new double[3];              //*dimensions of plane
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

//        short nSlices = 150; // Default to 50 slices
        // My interval is upper bounded by the greatest dimension since we need
        // to respect the interval for all dimensions (assumption we make).
        
        List<Integer> allDimensions = new LinkedList<Integer>();
        allDimensions.add(volume.getDimX());
        allDimensions.add(volume.getDimY());
        allDimensions.add(volume.getDimZ());
        Collections.sort(allDimensions);
        
        double interval = Math.sqrt( Math.pow(allDimensions.get(2), 2) + Math.pow(allDimensions.get(1), 2) ) / ControlOptions.N_SLICES;
        
        TFColor oldColor = TFColor.getWhiteBackground();
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                oldColor.toWhiteBackground();
                for (double t = interval*(ControlOptions.N_SLICES/-2); t <= interval*(ControlOptions.N_SLICES/2) ; t=t+interval ){
                    
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + viewVec[0] * t + volumeCenter[0];
                     pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + viewVec[1] * t + volumeCenter[1] ;
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + viewVec[2] * t + + volumeCenter[2] ;
                    
                    try {
                        int val = (int)InterpolationMethods.getVoxel(pixelCoord, volume);
                        TFColor transferColor = tFunc.getColor(val);
                        
                        voxelColor.r = oldColor.r + transferColor.a * ( transferColor.r - oldColor.r );
                        voxelColor.g = oldColor.g + transferColor.a * ( transferColor.g - oldColor.g );
                        voxelColor.b = oldColor.b + transferColor.a * ( transferColor.b - oldColor.b );

                        oldColor.r = voxelColor.r;
                        oldColor.g = voxelColor.g;
                        oldColor.b = voxelColor.b;

                    } catch (ArrayIndexOutOfBoundsException ex ){
                        /* If a voxel is out of bounds (meaning that we are going outside the limits, then we 
                        just ignore. This is possible because we take the maximum dimension as reference. */
                        
                    }                    
                }
                                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;     //multiply with 255 for the color
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor); //* you basically "store" the values in the finalimage that you are seeing
            }
        }

    }
    
    //-------------------------------------------------------------------------------------------------------------------
    void levoys(double[] viewMatrix) {

        // clear image
//        for (int j = 0; j < image.getHeight(); j++) {
//            for (int i = 0; i < image.getWidth(); i++) {
//                image.setRGB(i, j, 0);
//            }
//        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];           //*view vec has the values of the voxels for a ray (pen-> bottle 1,1,0 1,1,1 1,1,2),
        double[] uVec = new double[3];              //*dimensions of plane
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        TFColor voxelColor = new TFColor();
        
        List<Integer> allDimensions = new LinkedList<Integer>();
        allDimensions.add(volume.getDimX());
        allDimensions.add(volume.getDimY());
        allDimensions.add(volume.getDimZ());
        Collections.sort(allDimensions);
        
        double interval = Math.sqrt( Math.pow(allDimensions.get(2), 2) + Math.pow(allDimensions.get(1), 2) ) / ControlOptions.N_SLICES;
        double r = tfEditor2D.triangleWidget.radius;
        short fv = tfEditor2D.triangleWidget.baseIntensity;
        double minThreshold = tfEditor2D.triangleWidget.minThreshold;
        double maxThreshold = tfEditor2D.triangleWidget.maxThreshold;
        
        System.out.println("Using tmin: " + minThreshold);
        System.out.println("Using tmax: " + maxThreshold);
        
        TFColor transferColor = tfEditor2D.triangleWidget.color;
//        TFColor oldColor = TFColor.getBlackBackground();
        double [] shade = new double [3];
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                voxelColor.toBlackBackground();
                for (double t = interval*(ControlOptions.N_SLICES/-2); t <= interval*(ControlOptions.N_SLICES/2) ; t=t+interval ){
                    
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + viewVec[0] * t + volumeCenter[0];
                     pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + viewVec[1] * t + volumeCenter[1] ;
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + viewVec[2] * t + + volumeCenter[2] ;
                    
                    
                    try {
                        int fx = (int)InterpolationMethods.getVoxel(pixelCoord, volume);
                        
                        VoxelGradient gradInter = InterpolationMethods.getGradient(pixelCoord, gradients);
                        float gradMagInter = gradInter.mag;
                        double alpha;
                        if( gradMagInter == 0 && fx == fv && gradMagInter >= minThreshold && gradMagInter <= maxThreshold ){
                            alpha = 1;
                        } else if( gradMagInter > 0  && (fx  - r*gradMagInter) <= fv  && 
                                fv <=  (fx  + r*gradMagInter) 
                                && gradMagInter >= minThreshold && gradMagInter <= maxThreshold
                                ){
                            alpha = 1 - (1/r) * Math.abs( ((double)fv - (double)fx)/(double)gradMagInter );
                        }   
                        else{
                            alpha = 0;
                        }

                        if (ControlOptions.SHADDING) {
                            shade = getShade(viewMatrix, pixelCoord, transferColor, gradInter);
                        } else {
                            shade[0] = transferColor.r;
                            shade[1] = transferColor.g;
                            shade[2] = transferColor.b;
                        }


                        voxelColor.r = voxelColor.r + alpha * ( shade[0] - voxelColor.r );
                        voxelColor.g = voxelColor.g + alpha * ( shade[1] - voxelColor.g );
                        voxelColor.b = voxelColor.b + alpha * ( shade[2] - voxelColor.b );

//                        oldColor.r = voxelColor.r;
//                        oldColor.g = voxelColor.g;
//                        oldColor.b = voxelColor.b;

                    } catch (ArrayIndexOutOfBoundsException ex ){
                        /* If a voxel is out of bounds (meaning that we are going outside the limits, then we 
                        just ignore. This is possible because we take the maximum dimension as reference. */
                        
                    }                    
                }
                                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = 1 <= 1.0 ? (int) Math.floor(1 * 255) : 255;     //multiply with 255 for the color
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor); //* you basically "store" the values in the finalimage that you are seeing
            }
        }

    }

    public double [] getShade(double[] viewMatrix, double[] coord, TFColor iDiff, VoxelGradient gradient) {
        double[] N = new double[3];
        
//        VoxelGradient gradient = InterpolationMethods.getGradient(coord, gradients);
        VectorMath.setVector(N, gradient.x / gradient.mag, gradient.y / gradient.mag, gradient.z / gradient.mag);

        double[] V = new double[3];
        VectorMath.setVector(V, -viewMatrix[2], -viewMatrix[6], -viewMatrix[10]);
        double VLength = VectorMath.length(V);

        double[] L = new double[3];
        VectorMath.setVector(L, V[0] / VLength, V[1] / VLength, V[2] / VLength);
        
        double[] H = new double[3];
        VectorMath.setVector(H, 2*L[0], 2*L[1], 2*L[2]);

        double HLength = VectorMath.length(H);
        VectorMath.setVector(H, H[0] / HLength, H[1] / HLength, H[2] / HLength);

        double iAmb = 1; // white light
        double kAmb = 0.1;
//        double iDiff = 1;
        double kDiff = 0.7;
        double kSpec = 0.2;
        double alpha = 10;

        double shade = 0;
        double rShade = 0;
        double gShade = 0;
        double bShade = 0;

        double ambTerm = iAmb * kAmb;
        double c1 = VectorMath.dotproduct(N, L);
        double c2 = VectorMath.dotproduct(N, H);
//        if(c1 > 0 && c2 > 0){
//            shade = iAmb * kAmb + iDiff * kDiff * c1 + iDiff * kSpec * Math.pow(c2, alpha);
//        }
//        else {
//            shade = 1;
//        }
        if (ambTerm > 0) 
            shade += ambTerm;
        double rDiffTerm = iDiff.r * kDiff * c1;
        double gDiffTerm = iDiff.g * kDiff * c1;
        double bDiffTerm = iDiff.b * kDiff * c1;
//        double diffTerm = iDiff * kDiff * (VectorMath.dotproduct(N, L));
        if (rDiffTerm > 0)
            rShade += rDiffTerm;   
        if (gDiffTerm > 0)
            gShade += gDiffTerm; 
        if (bDiffTerm > 0)
            bShade += bDiffTerm; 
        double specTerm = kSpec * Math.pow(c2, alpha);
        if (specTerm > 0){
            rShade += specTerm;
            gShade += specTerm;
            bShade += specTerm;
        }
            
        double [] shades = {rShade, gShade, bShade};
        return shades;
    }

    
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        switch(ControlOptions.VISUALIZATION_MODE){
            case ControlOptions.SLICER_OPT:
                System.out.println("Using Slicer");
                this.slicer(viewMatrix);
            break;
            case ControlOptions.MIP_OPT:
                System.out.println("Using MIP");
                this.mip(viewMatrix);    
            break;
            case ControlOptions.COMPOSITING_OPT:
                System.out.println("Using Compositing");
                this.compositing(viewMatrix);
            break;
            case ControlOptions.TF2D_OPT:
                System.out.println("2D Transfer function");
                this.levoys(viewMatrix);
            break;
            default:
                this.slicer(viewMatrix);
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}