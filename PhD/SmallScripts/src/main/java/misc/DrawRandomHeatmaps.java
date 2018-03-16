/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import com.lowagie.text.DocumentException;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.graphics.Heatmap;

/**
 *
 * @author harmjan
 */
public class DrawRandomHeatmaps {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            Heatmap map = new Heatmap();

            int numSamples = 10;
            int numProbes = 10;
            int numTissues = 10;

            int width = numSamples * 500;
            int height = numSamples * 500;

            double[][] values = new double[numProbes][numSamples];
            double[][] values2 = new double[numProbes][numSamples];

            String[] rowheader = new String[numProbes];
            String[] rowheader2 = new String[numProbes];
            String[] colheader = new String[numSamples];

            for (int i = 0; i < numProbes; i++) {
                rowheader[i] = "Probe " + i;
                rowheader2[i] = "Cell type " + i;
                for (int j = 0; j < numSamples; j++) {
                    colheader[j] = "Sample " + j;
                    double d = Math.random();
                    values[i][j] = d;
                    
                    values2[i][j] = d*Math.random();
                }
            }
            String fileOut = "/Volumes/iSnackHD/tmp/plot.pdf";
            String fileOut2 = "/Volumes/iSnackHD/tmp/plot2.pdf";
            try {
                Heatmap.drawHeatmap(values, rowheader, colheader, width, height, fileOut, Heatmap.Output.PDF);
                Heatmap.drawHeatmap(values2, rowheader2, colheader, width, height, fileOut2, Heatmap.Output.PDF);
            } catch (DocumentException ex) {
                Logger.getLogger(DrawRandomHeatmaps.class.getName()).log(Level.SEVERE, null, ex);
            }

        } catch (IOException e) {
        }
    }
}
