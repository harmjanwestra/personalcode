/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.heterogeneity;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.math.stats.Heterogeneity;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CalculateHeterogeneity {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String[] datasetnames = new String[]{"EGCUT", "SHIP_TREND", "Groningen-HT12", "Groningen-H8v2", "Rotterdam", "DILGOM", "INCHIANTI", "HVH-HT12v3", "HVH-HT12v4"};
        // "EGCUT","SHIP_TREND","Groningen-HT12","Groningen-H8v2","Rotterdam","DILGOM","INCHIANTI","HVH-HT12v3","HVH-HT12v4";
        int[] datasetweights = new int[]{891, 963, 1240, 229, 762, 509, 611, 43, 63};
        // 891,963,1240,229,762,509,611,43,63


        try {
//            TextFile scriptOut = new TextFile("/Volumes/iSnackHD/SkyDrive/latesteQTLs/FilesForMETAL/script.mtl", TextFile.W);
//            scriptOut.writeln("VERBOSE ON");
//            scriptOut.writeln();



            TextFile in = new TextFile("/Volumes/iSnackHD/SkyDrive/latesteQTLs/transFDR0.05.txt.gz", TextFile.R);

            HashSet<String> visitedProbes = new HashSet<String>();

            TextFile out = new TextFile("/Volumes/iSnackHD/SkyDrive/latesteQTLs/transFDR0.05-WithHeteroGeneity.txt", TextFile.W);



            String[] data = in.readLineElems(TextFile.tab);
            out.writeln(Strings.concat(data, Strings.tab) + "\tISq\tISqPval");
            data = in.readLineElems(TextFile.tab);
            while (data != null) {
                if (!visitedProbes.contains(data[4])) {

                    Double metaZ = Double.parseDouble(data[eQTLTextFile.METAZ]);

                    String[] zscoreElems = data[eQTLTextFile.DATASETZSCORE].split(";");
                    if (zscoreElems.length == 1) {
                        zscoreElems = data[eQTLTextFile.DATASETZSCORE].split(",");
                    }
                    Double[] datasetZ = new Double[zscoreElems.length];

                    int nrDs = 0;
                    for (int d = 0; d < datasetZ.length; d++) {
                        String zScore = zscoreElems[d];
                        if (!zScore.equals("-")) {
                            // no p-value to calculate.. do not output pair
                            double z = Double.parseDouble(zScore);
                            datasetZ[d] = z;
                            nrDs++;
                        }
                    }
                    
                    

                    Pair<Double, Double> het = Heterogeneity.getISq(datasetZ, datasetweights);

                    out.writeln(Strings.concat(data, Strings.tab) + "\t" + het.getLeft() + "\t" + het.getRight());


                    // convert ZScore to pval..




//                        out.writeln(data[1] + "-" + data[4] + "\t" + alleleAssessed + "\t" + alleleElems[otherAllele] + "\t" + p + "\t" + z + "\t" + datasetweights[datasetNr]);
//                    }
                }
                visitedProbes.add(data[4]);
                data = in.readLineElems(TextFile.tab);



//                MARKER Pair
//ALLELE refallele nonrefallele
//EFFECT ZScore
//PVALUE PValue 
//WEIGHT Weight
//PROCESS inputfile1.txt

//                scriptOut.writeln(
//                        "MARKER Pair\n"
//                        + "ALLELE refallele nonrefallele\n"
//                        + "EFFECT ZScore\n"
//                        + "WEIGHT Weight\n"
//                        + "PROCESS /Volumes/iSnackHD/SkyDrive/latesteQTLs/FilesForMETAL/" + datasetnames[datasetNr] + ".tbl");
//                scriptOut.writeln("");
            }
//            scriptOut.writeln("OUTFILE METAANALYSIS.tbl");
//            scriptOut.writeln("ANALYZE HETEROGENEITY");
//            scriptOut.writeln("QUIT");
//            scriptOut.close();
// requires:
            // pair\trefallele nonrefallele\tpvalue\tzscore\tweight
            System.out.println("Done.");
            out.close();
            in.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
