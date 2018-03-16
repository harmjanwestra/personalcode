/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class AnandAndiappanCisEQTLFileCreator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        // create list of neutrophil specific cis-eQTLs
        // load the neutrophil specificity Z-score from the neutrophil specificity analysis
        // load the z-scores for these effects from the large eQTL meta
        String significanceFile = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-09-27-FDREstimates/CellTypeInteractionZScore-null.txt";
        String outfilename = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-11-19-AnandNeutrophils/RequestedEQTLs.txt";
        String eQTLfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/2012-12-21-CisAssociationsProbeLevelFDR0.5.txt";

        try {
            HashMap<String, Double> cellTypeSpecificityZScore = new HashMap<String, Double>();
            HashMap<String, Double> cellTypeSpecificityZScoreFDR = new HashMap<String, Double>();
            TextFile tf = new TextFile(significanceFile, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String snpprobe = elems[0];
                double ctz = Double.parseDouble(elems[2]);
                double fdr = Double.parseDouble(elems[4]);

                cellTypeSpecificityZScore.put(snpprobe, ctz);
                cellTypeSpecificityZScoreFDR.put(snpprobe, fdr);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile out = new TextFile(outfilename, TextFile.W);
            TextFile eqtl = new TextFile(eQTLfile, TextFile.R);
            out.writeln(eqtl.readLine() + "\tCellTypeSpecificityZScore\tCellTypeSpecificityZScoreFDR\tSpecificForCellType");
            elems = eqtl.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[1];
                String probe = elems[4];
                String snpprobe = snp + "-" + probe;

                Double ctz = cellTypeSpecificityZScore.get(snpprobe);
                if (ctz != null) {
                    Double ctzfdr = cellTypeSpecificityZScoreFDR.get(snpprobe);
                    String type = "Generic";
                    if (ctzfdr < 0.05) {
                        if (ctz > 0) {
                            type = "Neutrophil specific";
                        } else {
                            type = "Lymphocyte specific";
                        }
                    }
                    out.writeln(Strings.concat(elems, Strings.tab) + "\t" + ctz + "\t" + ctzfdr + "\t" + type);
                }

                elems = eqtl.readLineElems(TextFile.tab);
            }
            eqtl.close();
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
