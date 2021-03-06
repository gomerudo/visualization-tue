/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import implementation.InterpolationMethods;


/**
 *
 * @author michel
 */
public class GradientVolume {

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
    }

    public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }

    
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }

    public VoxelGradient getVoxel(int i) {
        return data[i];
    }

    public int getDimX() {
        return dimX;
    }

    public int getDimY() {
        return dimY;
    }

    public int getDimZ() {
        return dimZ;
    }

    private void compute() {

        float xG, yG, zG;
        float valueR, valueL;
        
        for( int i = 0; i < volume.getDimX(); i++ ){
            for( int j = 0; j < volume.getDimY(); j++ ){
                for( int k = 0; k < volume.getDimZ(); k++ ){
                    
                    /* For x */
                    valueL = i == 0 ? volume.getVoxel(i, j, k) : volume.getVoxel(i - 1, j, k);
                    valueR = i == volume.getDimX() - 1 ? volume.getVoxel(i, j, k) : volume.getVoxel(i + 1, j, k);
                    xG = (valueR - valueL)/2;
                    
                    /* For y */
                    valueL = j == 0 ? volume.getVoxel(i, j, k) : volume.getVoxel(i, j - 1, k);
                    valueR = j == volume.getDimY() - 1 ? volume.getVoxel(i, j, k) : volume.getVoxel(i, j + 1, k);
                    yG = (valueR - valueL)/2;
                    
                    /* For z */
                    valueL = k == 0 ? volume.getVoxel(i, j, k) : volume.getVoxel(i, j, k - 1);
                    valueR = k == volume.getDimZ() - 1 ? volume.getVoxel(i, j, k) : volume.getVoxel(i, j, k + 1);
                    zG = (valueR - valueL)/2;
                    
                    setGradient(i, j, k, new VoxelGradient(xG, yG, zG));
                }
            }
        }
        
    }
    
    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
}
