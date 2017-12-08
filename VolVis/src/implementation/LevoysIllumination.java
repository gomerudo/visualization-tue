/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package implementation;

import util.VectorMath;
import volume.VoxelGradient;
import volvis.TFColor;

/**
 *
 * @author gomerudo
 */
public class LevoysIllumination {
    
    public static double getAlpha(VoxelGradient gradInter, double minThreshold, double maxThreshold, int fx, short fv, double r){
        float gradMagInter = gradInter.mag;
        if( gradMagInter == 0 && fx == fv && // Normal condition 
                gradMagInter >= minThreshold && gradMagInter <= maxThreshold  // Extension of triangle widget
          ){
            return 1;
        } else if( gradMagInter > 0  && 
                (fx  - r*gradMagInter) <= fv  &&  fv <=  (fx  + r*gradMagInter) 
                && gradMagInter >= minThreshold && gradMagInter <= maxThreshold // Extension of triangle widget
                ){
            return 1 - (1/r) * Math.abs( ((double)fv - (double)fx)/(double)gradMagInter );
        }   
        else{
            return 0;
        }
    }
    
    public static double [] getShade(double[] V, double [] L, double [] H, double[] coord, TFColor iDiff, VoxelGradient gradient) {
        double[] N = new double[3];
        
        VectorMath.setVector(N, gradient.x / gradient.mag, gradient.y / gradient.mag, gradient.z / gradient.mag);

        double iAmb = 1; // white light
        double kAmb = 0.1;
//        double iDiff = 1;
        double kDiff = 0.7;
        double kSpec = 0.2;
        double alpha = 10;

        double rShade = 0;
        double gShade = 0;
        double bShade = 0;

        double ambTerm = iAmb * kAmb;
        double c1 = VectorMath.dotproduct(N, L);
        double c2 = VectorMath.dotproduct(N, H);

        if (ambTerm > 0) {
            gShade = bShade = rShade += ambTerm;
        }

        double rDiffTerm = iDiff.r * kDiff * c1;
        double gDiffTerm = iDiff.g * kDiff * c1;
        double bDiffTerm = iDiff.b * kDiff * c1;
        
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

    
}
