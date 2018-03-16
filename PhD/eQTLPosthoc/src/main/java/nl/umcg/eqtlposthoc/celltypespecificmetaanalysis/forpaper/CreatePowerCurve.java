/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.HashSet;
import java.util.LinkedList;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CreatePowerCurve {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        Meta m = new Meta();
        try {
            boolean flipEffects = true;
            boolean forcePresenceOfEQTL = true;
            String output = "/Volumes/iSnackHD/AeroFS/2013-12-06-IterativeMetaAllPermutations-EQTLsPresentInAllDatasets-Bonferroni/";
            Gpio.createDir(output);
            String probetranslationfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            String[] datasetFileDirs = new String[]{
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-29-DILGOM/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-08-19-INCHIANTI/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-08-19-KORA-F4/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-24_ROTTERDAM-STUDY/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-25-EGCUT/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-23_SHIP-TREND/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-12-GroningenHT12v3/",};

            // DILGOM, INCHIANTI, KORA, RS, EGCUT, SHIP, GRNG
            String[] datasetNames = new String[]{
                "DILGOM",
                "INCHIANTI",
                "KORA-F4",
                "Rotterdam",
                "EGCUT",
                "SHIP-TREND",
                "GroningenHT12v3"
            };

            int[] sampleSizes = new int[]{
                508,
                606,
                740,
                755,
                891,
                963,
                1220
            };

            boolean[] platformIsHT12v3 = new boolean[]{
                true,
                true,
                true,
                false,
                true,
                true,
                true
            };

            double zScoreThreshold = 2.605579924;// 4.445975 // bonferroni;

            HashSet<String> covariatesToTest = new HashSet<String>();
            covariatesToTest.add("CellTypeInteractionZScore");
            covariatesToTest.add("MainEffectZScore");

            TextFile tfOut = new TextFile(output + "Summary.txt", TextFile.W);
            tfOut.writeln("Iteration\tDatasets\tSampleSize\tNrGeneric\tNrLymphocyte\tNrNeutophil\tNrLymphocyteProbes\tNrNeutroProbes");
//            for (int iter = 0; iter < datasetNames.length; iter++) {
//                String[] datasetLocationsToInclude = new String[iter + 1];
//                String[] datasetNamesToInclude = new String[iter + 1];
//                boolean[] datasetPlatformisHT12v3ToInclude = new boolean[iter + 1];
//
//                for (int q = 0; q < iter + 1; q++) {
//                    datasetLocationsToInclude[q] = datasetFileDirs[q];
//                    datasetNamesToInclude[q] = datasetNames[q];
//                    datasetPlatformisHT12v3ToInclude[q] = platformIsHT12v3[q];
//                }
//
//                String outdir = output + "/Iteration" + (iter + 1) + "/";
//                m.run(outdir, probetranslationfile, datasetLocationsToInclude, datasetNamesToInclude, datasetPlatformisHT12v3ToInclude, false, flipEffects, covariatesToTest);
//
//                TextFile tf2 = new TextFile(outdir + "Vector-CellTypeInteractionZScore.txt", TextFile.R);
//                tf2.readLine();
//                String[] elems = tf2.readLineElems(TextFile.tab);
//                int nrLympho = 0;
//                int nrNeutro = 0;
//                int nrGeneric = 0;
//                tf2.readLine();
//                while (elems != null) {
//                    String eQTL = elems[0];
//                    Double d = Double.parseDouble(elems[1]);
//                    if (Math.abs(d) > zScoreThreshold) {
//                        if (d < 0) {
//                            nrLympho++;
//                        } else {
//                            nrNeutro++;
//                        }
//                    } else {
//                        nrGeneric++;
//                    }
//                    elems = tf2.readLineElems(TextFile.tab);
//                }
//                tf2.close();

//                tfOut.writeln((iter + 1) + "\t" + Strings.concat(datasetNamesToInclude, Strings.semicolon) + "\t" + 0 + "\t" + nrGeneric + "\t" + nrLympho + "\t" + nrNeutro);
//            }
            
            // randomize this stuff a bit.
            int permutation = 1;
            for (int numDatasets = 1; numDatasets < datasetNames.length + 1; numDatasets++) {
                LinkedList<String> list = Combinations.comb(numDatasets, datasetNames.length);
                for (String s : list) { // each permutation

                    String[] datasets = s.split(" ");

                    String[] datasetLocationsToInclude = new String[numDatasets];
                    String[] datasetNamesToInclude = new String[numDatasets];
                    boolean[] datasetPlatformisHT12v3ToInclude = new boolean[numDatasets];

                    int sampleSize = 0;
                    for (int iter = 0; iter < numDatasets; iter++) {
                        int dsId = Integer.parseInt(datasets[iter]);
                        datasetLocationsToInclude[iter] = datasetFileDirs[dsId];
                        datasetNamesToInclude[iter] = datasetNames[dsId];
                        datasetPlatformisHT12v3ToInclude[iter] = platformIsHT12v3[dsId];
                        sampleSize += sampleSizes[dsId];
                    }

                    String outdir = output + "/Iteration" + permutation + "/";
//                    m.run(outdir, probetranslationfile, datasetLocationsToInclude, datasetNamesToInclude, datasetPlatformisHT12v3ToInclude, false, flipEffects, covariatesToTest, forcePresenceOfEQTL);

                    TextFile tf2 = new TextFile(outdir + "Vector-CellTypeInteractionZScore.txt", TextFile.R);
                    tf2.readLine();
                    String[] elems = tf2.readLineElems(TextFile.tab);
                    int nrLympho = 0;
                    int nrNeutro = 0;
                    int nrGeneric = 0;
                    tf2.readLine();
                    HashSet<String> neutroprobes = new HashSet<String>();
                    HashSet<String> lymphoprobes = new HashSet<String>();
                    while (elems != null) {
                        String eQTL = elems[0];
                        String probe = eQTL.split("-")[1];
                        Double d = Double.parseDouble(elems[1]);
                        if (Math.abs(d) > zScoreThreshold) {
                            if (d < 0) {
                                nrLympho++;
                                lymphoprobes.add(probe);
                            } else {
                                nrNeutro++;
                                neutroprobes.add(probe);
                            }
                        } else {
                            nrGeneric++;
                        }
                        elems = tf2.readLineElems(TextFile.tab);
                    }
                    tfOut.writeln(permutation + "\t" + Strings.concat(datasetNamesToInclude, Strings.semicolon) + "\t" + sampleSize + "\t" + nrGeneric + "\t" + nrLympho + "\t" + nrNeutro+ "\t" + lymphoprobes.size() + "\t" + neutroprobes.size());
                    tf2.close();
                    permutation++;
                }
            }
            tfOut.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
