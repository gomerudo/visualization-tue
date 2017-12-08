/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package implementation;

import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author gomerudo
 */
public class InterpolationMethods {
    
    public static double getVoxel(double[] coord, Volume volume) {          //interpolated values for MIP
        
        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()  //same as slicers gett value
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }
        
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        
        
        if ( x > volume.getDimX()          //with this i dont get out of bounds error
                || y > volume.getDimY()
                || z > volume.getDimZ() ) {
            return 0;
        }
        
        int x0 = (int) Math.floor(x);
        int y0 = (int) Math.floor(y);
        int z0 = (int) Math.floor(z);  
        
     
        int x1 = (int) Math.ceil(x);
        int y1 = (int) Math.ceil(y);
        int z1 = (int) Math.ceil(z); 
        
        

        if (x1 < 0 || x1 > volume.getDimX() || y1 < 0 || y1 > volume.getDimY()
                || z1 < 0 || z1 > volume.getDimZ()) {
            return 0;
        }
        
        double alpha = x - x0;
        double beta = y - y0;
        double gama = z - z0;
        
        double A = volume.getVoxel(x0, y0, z0);
        double B = volume.getVoxel(x1, y0, z0);
        double C = volume.getVoxel(x0, y1, z0);
        double D = volume.getVoxel(x1, y1, z0);
        double E = volume.getVoxel(x0, y0, z1);
        double F = volume.getVoxel(x1, y0, z1);
        double G = volume.getVoxel(x0, y1, z1);
        double H = volume.getVoxel(x1, y1, z1);
        
        double interVal = (double) (1 - gama)*( (1 - alpha)*(A + beta*(C - A)) + alpha*( beta*(D - B ) + B)  ) 
                + beta*( gama*( alpha*(H - G) + G ) )
                + (1 - beta)*( gama*( E + alpha*(F - E) ) )
                ;
//        double interVal = (double) ((1 - alpha) * (1 - beta) * (1 - gama) * volume.getVoxel(x0, y0, z0)
//                + alpha * (1 - beta) * (1 - gama) * volume.getVoxel(x1, y0, z0)
//                + (1 - alpha) * beta * (1 - gama) * volume.getVoxel(x0, y1, z0)
//                + alpha * beta * (1 - gama) * volume.getVoxel(x1, y1, z0)
//                + (1 - alpha) * (1 - beta) * gama * volume.getVoxel(x0, y0, z1)
//                + alpha * (1 - beta) * gama * volume.getVoxel(x1, y0, z1)
//                + (1 - alpha) * beta * gama * volume.getVoxel(x0, y1, z1)
//                + alpha * beta * gama * volume.getVoxel(x1, y1, z1));

        return interVal;
          
    }

    public static VoxelGradient getGradient(double[] coord, GradientVolume volume) {          //interpolated values for MIP
        
        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()  //same as slicers gett value
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return new VoxelGradient();
        }
        
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        
        
        if ( x > volume.getDimX()          //with this i dont get out of bounds error
                || y > volume.getDimY()
                || z > volume.getDimZ() ) {
            return new VoxelGradient();
        }
        
        int x0 = (int) Math.floor(x);
        int y0 = (int) Math.floor(y);
        int z0 = (int) Math.floor(z);  
        
     
        int x1 = (int) Math.ceil(x);
        int y1 = (int) Math.ceil(y);
        int z1 = (int) Math.ceil(z); 
        
        

        if (x1 < 0 || x1 > volume.getDimX() || y1 < 0 || y1 > volume.getDimY()
                || z1 < 0 || z1 > volume.getDimZ()) {
            return new VoxelGradient();
        }
        
        double alpha = x - x0;
        double beta = y - y0;
        double gama = z - z0;
        
        
        double A = volume.getGradient(x0, y0, z0).x;
        double B = volume.getGradient(x1, y0, z0).x;
        double C = volume.getGradient(x0, y1, z0).x;
        double D = volume.getGradient(x1, y1, z0).x;
        double E = volume.getGradient(x0, y0, z1).x;
        double F = volume.getGradient(x1, y0, z1).x;
        double G = volume.getGradient(x0, y1, z1).x;
        double H = volume.getGradient(x1, y1, z1).x;
        
        double interValX = (double) (1 - gama)*( (1 - alpha)*(A + beta*(C - A)) + alpha*( beta*(D - B ) + B)  ) 
                + beta*( gama*( alpha*(H - G) + G ) )
                + (1 - beta)*( gama*( E + alpha*(F - E) ) )
                ;

        
//        double interValX = (double) ((1 - alpha) * (1 - beta) * (1 - gama) * volume.getGradient(x0, y0, z0).x
//                + alpha * (1 - beta) * (1 - gama) * volume.getGradient(x1, y0, z0).x
//                + (1 - alpha) * beta * (1 - gama) * volume.getGradient(x0, y1, z0).x
//                + alpha * beta * (1 - gama) * volume.getGradient(x1, y1, z0).x
//                + (1 - alpha) * (1 - beta) * gama * volume.getGradient(x0, y0, z1).x
//                + alpha * (1 - beta) * gama * volume.getGradient(x1, y0, z1).x
//                + (1 - alpha) * beta * gama * volume.getGradient(x0, y1, z1).x
//                + alpha * beta * gama * volume.getGradient(x1, y1, z1).x
//                );
        
        A = volume.getGradient(x0, y0, z0).y;
        B = volume.getGradient(x1, y0, z0).y;
        C = volume.getGradient(x0, y1, z0).y;
        D = volume.getGradient(x1, y1, z0).y;
        E = volume.getGradient(x0, y0, z1).y;
        F = volume.getGradient(x1, y0, z1).y;
        G = volume.getGradient(x0, y1, z1).y;
        H = volume.getGradient(x1, y1, z1).y;

        double interValY = (double) (1 - gama)*( (1 - alpha)*(A + beta*(C - A)) + alpha*( beta*(D - B ) + B)  ) 
                + beta*( gama*( alpha*(H - G) + G ) )
                + (1 - beta)*( gama*( E + alpha*(F - E) ) )
                ;
//        double interValY = (double) ((1 - alpha) * (1 - beta) * (1 - gama) * volume.getGradient(x0, y0, z0).y
//                + alpha * (1 - beta) * (1 - gama) * volume.getGradient(x1, y0, z0).y
//                + (1 - alpha) * beta * (1 - gama) * volume.getGradient(x0, y1, z0).y
//                + alpha * beta * (1 - gama) * volume.getGradient(x1, y1, z0).y
//                + (1 - alpha) * (1 - beta) * gama * volume.getGradient(x0, y0, z1).y
//                + alpha * (1 - beta) * gama * volume.getGradient(x1, y0, z1).y
//                + (1 - alpha) * beta * gama * volume.getGradient(x0, y1, z1).y
//                + alpha * beta * gama * volume.getGradient(x1, y1, z1).y
//                );
        
        A = volume.getGradient(x0, y0, z0).z;
        B = volume.getGradient(x1, y0, z0).z;
        C = volume.getGradient(x0, y1, z0).z;
        D = volume.getGradient(x1, y1, z0).z;
        E = volume.getGradient(x0, y0, z1).z;
        F = volume.getGradient(x1, y0, z1).z;
        G = volume.getGradient(x0, y1, z1).z;
        H = volume.getGradient(x1, y1, z1).z;

        double interValZ = (double) (1 - gama)*( (1 - alpha)*(A + beta*(C - A)) + alpha*( beta*(D - B ) + B)  ) 
                + beta*( gama*( alpha*(H - G) + G ) )
                + (1 - beta)*( gama*( E + alpha*(F - E) ) )
                ;

//        double interValZ = (double) ((1 - alpha) * (1 - beta) * (1 - gama) * volume.getGradient(x0, y0, z0).z
//                + alpha * (1 - beta) * (1 - gama) * volume.getGradient(x1, y0, z0).z
//                + (1 - alpha) * beta * (1 - gama) * volume.getGradient(x0, y1, z0).z
//                + alpha * beta * (1 - gama) * volume.getGradient(x1, y1, z0).z
//                + (1 - alpha) * (1 - beta) * gama * volume.getGradient(x0, y0, z1).z
//                + alpha * (1 - beta) * gama * volume.getGradient(x1, y0, z1).z
//                + (1 - alpha) * beta * gama * volume.getGradient(x0, y1, z1).z
//                + alpha * beta * gama * volume.getGradient(x1, y1, z1).z
//                );

        return new VoxelGradient((float)interValX, (float)interValY, (float)interValZ);
          
    }

}
