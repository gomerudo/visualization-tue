/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package util;

/**
 *
 * @author gomerudo
 */
public class ControlOptions {
    public static int VISUALIZATION_MODE = 1;
    
    public static final int SLICER_OPT = 1;
    public static final int MIP_OPT = 2;
    public static final int COMPOSITING_OPT = 3;
    public static final int TF2D_OPT = 4;
    
    public static boolean SHADDING = false;
    public static boolean LOW_RESOLUTION = false; // Not yet implemented
    
    public static int N_SLICES = 50;
    public static int BKP_SLICES = 50;
    
    public static double I_AMB = 1;
    public static double K_AMB = 0.1;
    public static double K_DIFF = 0.7;
    public static double K_SPEC = 0.2;
    public static double SHADE_ALPHA = 10;    
    
}
