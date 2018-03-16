/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class DetermineNumberOfDatasetsPerSNPProbePair {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here

            DetermineNumberOfDatasetsPerSNPProbePair.run("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PostQC/eQTLProbesFDR0.05.txt");
        } catch (IOException ex) {
            Logger.getLogger(DetermineNumberOfDatasetsPerSNPProbePair.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static void run(String file) throws IOException {
        TextFile tf = new TextFile(file, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        String[] elems = tf.readLineElems(TextFile.tab);



        HashMap<Integer, Integer> oppositeSampleCtr = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> sameDirectionSampleCtr = new HashMap<Integer, Integer>();

        HashMap<Pair<Integer, Integer>, Integer> dsSameDirectionCounter = new HashMap<Pair<Integer, Integer>, Integer>();
        for (int i = 1; i < 10; i++) {
            for (int j = 0; j < i + 1; j++) {
                dsSameDirectionCounter.put(new Pair<Integer, Integer>(i, j), 0);
            }
        }

        int[] dsCounter = new int[15];
        int[] oppositeDsCounter = new int[15];
        int[] sameDirectionDsCounter = new int[15];


        TextFile out = new TextFile(file + "-NrSamplesWithIdenticalEffectDirection.txt", TextFile.W);
        int nrPositiveEQTLs = 0;
        int nrNegativeEQTLs = 0;
        while (elems != null) {

            String zscores = elems[eQTLTextFile.DATASETZSCORE];
            String sampleSizeStr = elems[13];
            String[] zscoreElems = zscores.split(";");
            if (zscoreElems.length == 1) {
                zscoreElems = zscores.split(",");
            }

            String[] sampleSizes = sampleSizeStr.split(";");
            if (sampleSizes.length == 1) {
                sampleSizes = sampleSizeStr.split(",");
            }

            int dsctr = 0;

            int dsWithSameDirection = 0;
            int dsWithOppositeSameDirection = 0;
            int nrSamplesWithSameDirection = 0;
            int nrSamplesWithOppositeDirection = 0;

            int nrPositive = 0;
            int nrNegative = 0;

            for (int i = 0; i < zscoreElems.length; i++) {
                if (!zscoreElems[i].equals("-")) {
                    dsctr++;
                    double z = Double.parseDouble(zscoreElems[i]);
                    int sampleSize = Integer.parseInt(sampleSizes[i]);
                    if (z >= 0) {
                        nrPositive += sampleSize;
                    } else {
                        nrNegative += sampleSize;
                    }
                }
            }

            boolean directionIsPositive = false;
            if (nrPositive >= nrNegative) {
                directionIsPositive = true;
            }

            for (int i = 0; i < zscoreElems.length; i++) {
                if (!zscoreElems[i].equals("-")) {
                    double z = Double.parseDouble(zscoreElems[i]);
                    int sampleSize = Integer.parseInt(sampleSizes[i]);
                    if ((z >= 0 && directionIsPositive) || (z < 0 && !directionIsPositive)) {
                        dsWithSameDirection++;
                        nrSamplesWithSameDirection += sampleSize;
                    } else {
                        dsWithOppositeSameDirection++;
                        nrSamplesWithOppositeDirection += sampleSize;
                    }
                }
            }


            if(directionIsPositive){
                nrPositiveEQTLs++;
            } else {
                nrNegativeEQTLs++;
            }

            oppositeDsCounter[dsWithOppositeSameDirection]++;
            sameDirectionDsCounter[dsWithSameDirection]++;
            dsCounter[dsctr]++;

            Pair<Integer, Integer> p = new Pair<Integer, Integer>(dsctr, dsWithSameDirection);
            int c = dsSameDirectionCounter.get(p);
            c++;
            dsSameDirectionCounter.put(p, c);

            Integer tmpCtr = oppositeSampleCtr.get(nrSamplesWithOppositeDirection);
            if (tmpCtr == null) {
                tmpCtr = 1;
            } else {
                tmpCtr++;
            }
            oppositeSampleCtr.put(nrSamplesWithOppositeDirection, tmpCtr);

//            if (dsWithSameDirection == 1) {
            System.out.println(directionIsPositive + "\t" + nrSamplesWithSameDirection + "\t" + nrSamplesWithOppositeDirection + "\t" + Strings.concat(elems, Strings.tab));
//            }

            out.writeln(directionIsPositive + "\t" + nrSamplesWithSameDirection + "\t" + nrSamplesWithOppositeDirection + "\t" + Strings.concat(elems, Strings.tab));
            tmpCtr = sameDirectionSampleCtr.get(nrSamplesWithSameDirection);
            if (tmpCtr == null) {
                tmpCtr = 1;
            } else {
                tmpCtr++;
            }
            sameDirectionSampleCtr.put(nrSamplesWithSameDirection, tmpCtr);



            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        out.close();

        System.out.println("NrPos: "+nrPositiveEQTLs);
        System.out.println("NrNeg: "+nrNegativeEQTLs);
        System.out.println("");
        System.out.println("Stats");
        System.out.println("-----");
        System.out.println("NrDatasets\tNrEQTLsWData\tNrDsSameDirection\tNrDSOppositeDirection");
        for (int i = 0; i < 15; i++) {
            System.out.println(i + "\t" + dsCounter[i] + "\t" + sameDirectionDsCounter[i] + "\t" + oppositeDsCounter[i]);
        }
//        System.out.println("");
//        System.out.println("NrSamplesVsEQTLsWithSameDirection");
//        System.out.println("---------------------------------");
//        Set<Entry<Integer,Integer>> entries = sameDirectionSampleCtr.entrySet();
//        for(Entry<Integer, Integer> e: entries){
//            System.out.println(e.getKey()+"\t"+e.getValue());
//        }
//        
        System.out.println("");
//        System.out.println("NrSamplesVsEQTLsWithOppositeDirection");
//        System.out.println("---------------------------------");
//        entries = oppositeSampleCtr.entrySet();
//        for(Entry<Integer, Integer> e: entries){
//            System.out.println(e.getKey()+"\t"+e.getValue());
//        }

        for (int i = 1; i < 10; i++) {
            for (int j = 0; j < i + 1; j++) {
                System.out.println(
                        i + "\t"
                        + j + "\t"
                        + (i - j) + "\t"
                        + dsSameDirectionCounter.get(new Pair<Integer, Integer>(i, j)));
            }
        }
    }
}
