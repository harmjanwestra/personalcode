/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import eqtlmappingpipeline.binarymeta.meta.graphics.ZScorePlot;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author harm-jan
 */
public class PlotZScores {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        // trans
        String[] dsNames = new String[]{"EGCUT", "SHIP-TREND", "Fehrmann-HT12", "Fehrmann-H8v2", "Rotterdam-Study", "DILGOM", "INCHIANTI", "HVH-HT12v3", "HVH-HT12v4", "Meta-Analysis"};
        String efilename = "D:\\SkyDrive\\latesteQTLs\\transFDR0.05.txt";
        String outname = "D:\\SkyDrive\\latesteQTLs\\transFDR0.05-ZScoreComparison.png";
        String sep = ";";
        // cis
//        String[] dsNames = new String[]{"EGCUT", "SHIP-TREND", "Fehrmann-HT12", "Fehrmann-H8v2", "Rotterdam-Study", "DILGOM", "INCHIANTI", "HVH-HT12v3", "HVH-HT12v4", "Meta-Analysis"};
//        String efilename = "D:\\SkyDrive\\latesteQTLs\\cisProbeLevelFDR0.05.txt";
//        String outname = "D:\\SkyDrive\\latesteQTLs\\cisProbeLevelFDR0.05-ZScoreComparison.png";
//        String sep = ",";
        PlotZScores p = new PlotZScores();
        try {
            p.plot(efilename, dsNames, outname, sep);
        } catch (IOException ex) {
            Logger.getLogger(PlotZScores.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(PlotZScores.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void plot(String eQTLFileName, String[] dsNames, String outfileName, String sep) throws IOException, Exception {
        ZScorePlot zsp = new ZScorePlot();
        zsp.init(dsNames.length, dsNames, false, outfileName);

        TextFile tf = new TextFile(eQTLFileName, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            Double metaZ = Double.parseDouble(elems[eQTLTextFile.METAZ]);
            String dsZscores = elems[eQTLTextFile.DATASETZSCORE];

            String[] dsZScoreArr = dsZscores.split(sep);

            for (int i = 0; i < dsZScoreArr.length; i++) {
                try {
                    Double z1 = Double.parseDouble(dsZScoreArr[i]);
                    for (int j = i + 1; j < dsZScoreArr.length; j++) {
                        try {
                            Double z2 = Double.parseDouble(dsZScoreArr[j]);
                            zsp.draw(z1, z2, i, j);
                        } catch (NumberFormatException e) {
                            
                        }
                    }
                    zsp.draw(z1, metaZ, i, dsNames.length - 1);
                } catch (NumberFormatException e) {
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }


        tf.close();

        zsp.write(outfileName);

    }
}
