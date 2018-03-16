/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import umcg.genetica.graphics.ScatterPlot;

/**
 *
 * @author harmjan
 */
public class ScatterplotTester {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here


        double[] x = new double[]{-4,-1,4,5,7};
        double[] y = new double[]{4,-1,-4,-5,-1100};
        ScatterPlot plot = new ScatterPlot(300, 300, x, y, ScatterPlot.OUTPUTFORMAT.PDF, "/Volumes/iSnackHD/tmp/test.png");

    }
}
