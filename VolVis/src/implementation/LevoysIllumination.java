/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package implementation;

import util.ControlOptions;
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
    
    public static double [] getShade(double [] L, double [] H, double[] coord, TFColor iDiff, VoxelGradient gradient) {
        double[] N = new double[3];
        
        VectorMath.setVector(N, gradient.x / gradient.mag, gradient.y / gradient.mag, gradient.z / gradient.mag);

        double rShade = 0;
        double gShade = 0;
        double bShade = 0;

        double ambTerm = ControlOptions.I_AMB * ControlOptions.K_AMB;
        double c1 = VectorMath.dotproduct(N, L);
        double c2 = VectorMath.dotproduct(N, H);

        double rDiffTerm = 0;
        double gDiffTerm = 0;
        double bDiffTerm = 0;
        
        // Formulas only apply when dot products are possitive
        if( c1 < 0 || c2 < 0){
            return new double []{iDiff.r, iDiff.g, iDiff.b};
        }

        // Just do computations when worth it (save runtime as much as possible)
        if(c1 > 0 || c2 > 0){
            double kDc1 = ControlOptions.K_DIFF * c1;
            rDiffTerm = iDiff.r * kDc1;
            gDiffTerm = iDiff.g * kDc1;
            bDiffTerm = iDiff.b * kDc1;
        }
        
        if (ambTerm > 0) {
            gShade = bShade = rShade += ambTerm; // All values start with the same factor
        }

        
        if (rDiffTerm > 0)
            rShade += rDiffTerm;   
        if (gDiffTerm > 0)
            gShade += gDiffTerm; 
        if (bDiffTerm > 0)
            bShade += bDiffTerm; 
        
        double specTerm = ControlOptions.K_SPEC * Math.pow(c2, ControlOptions.SHADE_ALPHA);
        if (specTerm > 0){
            rShade += specTerm;
            gShade += specTerm;
            bShade += specTerm;
        }
            
        double [] shades = {rShade, gShade, bShade};
        return shades;
    }

    
}
