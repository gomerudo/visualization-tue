/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

/**
 *
 * @author michel
 */
public class TFColor {
    public double r, g, b, a;

    public TFColor() {
        r = g = b = a = 1.0;
    }
    
    public TFColor(double red, double green, double blue, double alpha) {
        r = red;
        g = green;
        b = blue;
        a = alpha;
    }
    
    public void toWhiteBackground(){
        r = g = b = 1;
        a = 1;
    }
    
    public void toBlackBackground(){
        r = g = b = 0;
        a = 1;
    }
    
    public static TFColor getWhiteBackground(){
        return new TFColor(1, 1, 1, 1);
    }
    
    public static TFColor getBlackBackground(){
        return new TFColor(0, 0, 0, 1);
    }
            
    public boolean isWhiteBackground(){
        return this.r == 1 && this.g == 1 && this.b == 1;
    }
    
    public boolean isBlackBackground(){
        return this.r == 1 && this.g == 1 && this.b == 1;
    }
    @Override
    public String toString() {
        String text = "(" + r + ", " + g + ", " + b + ", " + a + ")";
        return text;
    }
}
