/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CreateReplicationTable {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        String metaAnalysisTable = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/2013-11-27-MetaAnalysisZScorePlots/ComparisonTable.txt";
        String replicationTable = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-11-27-ReplicationInCellTypeSpecificDatasets/CellTypeInteractionZScore-Summary.txt";
        String outfile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/2013-11-ManuscriptNatureMethods-Draft1/Supplementary/Supplementary Table 2 - Replication Effect sizes and mean expression.txt";
        String moreAnnotation = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/2013-11-ManuscriptNatureMethods-Draft1/Supplementary/Supplementary Table 1 - Results of cell type specificity analysis.txt";

        try {
            TextFile tf = new TextFile(replicationTable, TextFile.R);
            String[] header = tf.readLineElems(TextFile.tab);
            int nrdatasets = (header.length - 1) / 2;
            int nrEQTLs = tf.countLines();
            double[][] eQTLReplicationCorrelationCoefficients = new double[nrEQTLs][nrdatasets];
            tf.close();
            tf.open();
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            int ctr = 0;
            HashMap<String, Integer> eQTLToId = new HashMap<String, Integer>();

            while (elems != null) {
                String eQTL = elems[0];
                for (int d = 0; d < nrdatasets; d++) {
                    double z = Double.parseDouble(elems[(d * 2) + 1]);
                    eQTLReplicationCorrelationCoefficients[ctr][d] = z;
                }
                eQTLToId.put(eQTL, ctr);
                ctr++;
                elems = tf.readLineElems(TextFile.tab);

            }
            tf.close();

            HashMap<String, String> eQTLToGene = new HashMap<String, String>();
            HashMap<String, String> eQTLToPValue = new HashMap<String, String>();
            HashMap<String, String> eQTLToClassification = new HashMap<String, String>();
            TextFile tf3 = new TextFile(moreAnnotation, TextFile.R);
            tf3.readLine();
            String[] tf3elems = tf3.readLineElems(TextFile.tab);
            while (tf3elems != null) {
                String eQTL = tf3elems[1] + "-" + tf3elems[4];
                String gene = tf3elems[tf3elems.length - 5];
                String classification = tf3elems[tf3elems.length - 1];

                eQTLToPValue.put(eQTL, tf3elems[tf3elems.length - 2]);
                eQTLToGene.put(eQTL, gene);
                eQTLToClassification.put(eQTL, classification);
                tf3elems = tf3.readLineElems(TextFile.tab);
            }
            tf3.close();

            TextFile tfout = new TextFile(outfile, TextFile.W);
            String outheader = "SNP\tProbe\tGene\tClassification\tMetaAnalysisInteraction-R";
            for (int i = 0; i < nrdatasets; i++) {
                outheader += "\t" + header[1 + (i * 2)];
            }

            outheader += "\tAverageLymphoid\tAverageMyeloid";
            tfout.writeln(outheader);
            TextFile tf2 = new TextFile(metaAnalysisTable, TextFile.R);
            tf2.readLine(); // skip the headererer
            String[] tf2elems = tf2.readLineElems(TextFile.tab);

            HashMap<String, String> eQTLtoR = new HashMap<String, String>();

            while (tf2elems != null) {

                String eQTL = tf2elems[0];
                String metaZ = tf2elems[2];
                eQTLtoR.put(eQTL, metaZ);
                Integer eQTLId = eQTLToId.get(eQTL);
                String[] eQTLElems = eQTL.split("-");
                if (eQTLElems.length == 1) {
                    System.err.println("ERROR could not split eQTL: " + eQTL);
                } else {
                    String gene = eQTLToGene.get(eQTL);
                    if (gene == null) {
                        gene = "-";
                    }

                    String classification = eQTLToClassification.get(eQTL);

                    String outputln = eQTLElems[0] + "\t" + eQTLElems[1] + "\t" + gene + "\t" + classification + "\t" + ZScores.zToP(Double.parseDouble(metaZ));

                    double myeloidsum = 0;
                    double lymphoidsum = 0;
                    int lymphoidcounter = 0;
                    int myeloidcounter = 0;
                    if (eQTLId == null) {
                        System.err.println("ERROR: eQTL not found in replication: " + eQTL);
                    } else {

                        for (int d = 0; d < eQTLReplicationCorrelationCoefficients[eQTLId].length; d++) {
                            if (Double.isNaN(eQTLReplicationCorrelationCoefficients[eQTLId][d])) {
                                outputln += "\t-";
                            } else {
                                if (d < 4) {
                                    lymphoidsum += Math.abs(eQTLReplicationCorrelationCoefficients[eQTLId][d]);
                                    lymphoidcounter++;
                                } else {
                                    myeloidcounter++;
                                    myeloidsum += Math.abs(eQTLReplicationCorrelationCoefficients[eQTLId][d]);
                                }
                                outputln += "\t" + Math.abs(eQTLReplicationCorrelationCoefficients[eQTLId][d]);
                            }
                        }
                    }
                    if (lymphoidcounter > 0) {
                        lymphoidsum /= lymphoidcounter;
                    } else {
                        lymphoidsum = Double.NaN;
                    }
                    if (myeloidcounter > 0) {
                        myeloidsum /= myeloidcounter;
                    } else {
                        myeloidsum = Double.NaN;
                    }
                    tfout.writeln(outputln + "\t" + lymphoidsum + "\t" + myeloidsum);
                }
                tf2elems = tf2.readLineElems(TextFile.tab);
            }

            tf2.close();
            tfout.close();

            tf3.open();
            TextFile tf4 = new TextFile(tf3.getFileName() + "-WithR.txt", TextFile.W);
            tf4.writeln(tf3.readLine() + "\tInteractionPValue");
            tf3elems = tf3.readLineElems(TextFile.tab);
            while (tf3elems != null) {
                String eqtl = tf3elems[1] + "-" + tf3elems[4];
                tf4.writeln(Strings.concat(tf3elems, Strings.tab) + "\t" + ZScores.zToP(Double.parseDouble(eQTLtoR.get(eqtl))));
                tf3elems = tf3.readLineElems(TextFile.tab);
            }
            tf3.close();
            tf4.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
