/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import umcg.genetica.graphics.ViolinBoxPlot;
import umcg.genetica.math.stats.WilcoxonMannWhitney;

/**
 *
 * @author harmjan
 */
public class ViolinBoxPlotTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {
            double[][][] data = new double[0][0][0]; // format: vals[dataset][category][value] in your case: double[][][] data = vals[1][2][0], for 1 dataset, 2 categories, and x values (jagged array). you can just put the gene expression values in there
            String[] datasetnames = new String[0]; // the dataset names. in your case "GEO DATA" or something
            String[][] xLabels = new String[0][0]; // format: xlabels2[dataset][category] in your case xlabels[0][1] = "neutrophils" and xlabels[0][1] = "macrophages" or vice versa
            String ylabel = "relative expression";
            String outputfilename = "output.pdf";
            ViolinBoxPlot vbp = new ViolinBoxPlot();
            vbp.draw(data, datasetnames, xLabels, ylabel, ViolinBoxPlot.Output.PDF, outputfilename);
            
        } catch (Exception e) {
            
        }
        
        WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();
        // TODO code application logic here
//        try {
//            double[][][] data = new double[3][4][1000];
//            String[] datsetnames = new String[3];
//            String[][] xLabels1 = new String[3][4];
//            for (int x = 0; x < data.length; x++) {
//                datsetnames[x] = "Dataset-" + x;
//                double xfactor = 1;
//                if(x == 1){
//                    xfactor = -1;
//                }
//                for (int y = 0; y < data[x].length; y++) {
//                    xLabels1[x][y] = "Data-" + y;
//                    for (int z = 0; z < data[x][y].length; z++) {
//                        data[x][y][z] = (Math.random());
//                    }
//                }
//            }
//
            
//            vbp.draw(data, datsetnames, xLabels1, ViolinBoxPlot.Output.PDF, "/Volumes/iSnackHD/tmp/plottest.pdf");
//
//            
//
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
    }
}
