/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.endophenotypes;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CombineEndoPhenotypeCorrelationsWithTransEQTLTable {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
//         1980598TODO code application logic here
        try {
            int[] sampleSizes = new int[]{891, 963, 1240, 229, 762, 509, 611, 43, 63};
            String transFile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/PhD/2011-2012-eQTLMeta/Pack-Final/Supp/Supplementary Table 2 - Significant trans-eQTLs (FDR0.5).txt";
            String endofile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypes/Meta/MetaMerged-40PCsTransEQTLs.txt";
            String egcutcorrected = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypes/2013-02-EGCUTTransCorrectedForCovariates/Corrected/eQTLs.txt.gz";
            String egcutuncorrected = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypes/2013-02-EGCUTTransCorrectedForCovariates/UnCorrected/eQTLs.txt.gz";
            String outfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/PhD/2011-2012-eQTLMeta/Pack-Final/Supp/Supplementary Table 2 - Significant trans-eQTLs (FDR0.5)-WithEndoPhenotypeR2AndEGCUTRSquares.txt";
            String probeTranslationFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            CombineEndoPhenotypeCorrelationsWithTransEQTLTable c = new CombineEndoPhenotypeCorrelationsWithTransEQTLTable();
            c.run(endofile, transFile, probeTranslationFile, sampleSizes, outfile, egcutcorrected, egcutuncorrected);
        } catch (IOException e) {
            e.printStackTrace();
        }
//        CombineEndoPhenotypeCorrelationsWithTransEQTLTable c = new CombineEndoPhenotypeCorrelationsWithTransEQTLTable();
//        c.test();
    }

    public void run(String endoFile, String transEQTLFile, String probeTranslationFile, int[] sampleSizes, String outfile, String corrected, String uncorrected) throws IOException {

        TextFile tf1 = new TextFile(corrected, TextFile.R);
        TextFile tf2 = new TextFile(uncorrected, TextFile.R);
        HashMap<String, String> correctedZSCore = new HashMap<String, String>();
        HashMap<String, String> uncorrectedZSCore = new HashMap<String, String>();

        tf1.readLine();
        tf2.readLine();

        String[] elems = tf1.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            String probe = elems[4];
            String z = elems[eQTLTextFile.DATASECORR];
            double r = Double.parseDouble(z);
            z = "" + (r * r);
            correctedZSCore.put(snp + "-" + probe, z);
            elems = tf1.readLineElems(TextFile.tab);
        }

        System.out.println(correctedZSCore.size() + " corrected ZScores");


        elems = tf2.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            String probe = elems[4];
            String z = elems[eQTLTextFile.DATASECORR];
            double r = Double.parseDouble(z);
            z = "" + (r * r);
            uncorrectedZSCore.put(snp + "-" + probe, z);
            elems = tf2.readLineElems(TextFile.tab);
        }

        System.out.println(uncorrectedZSCore.size() + " corrected ZScores");
        tf1.close();
        tf2.close();


        ArrayList<String> toAdd = new ArrayList<String>();
        HashMap<String, Integer> probeToLn = new HashMap<String, Integer>();

        TextFile tf = new TextFile(endoFile, TextFile.R);

        elems = tf.readLineElems(TextFile.tab); // skip header
        String newHeader = "";
        for (int i = 6; i < elems.length; i += 6) {
            String headerElem = elems[i + 4];
            headerElem = headerElem.replace("Z-40PCs-", "");
            headerElem += " R2";
            newHeader += "\t" + headerElem;
        }
        elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;
        while (elems != null) {
            String probe = elems[0];
            String ln = "";
//            System.out.println(row);
            for (int i = 6; i < elems.length; i += 6) {
//                double origP = Double.parseDouble(elems[i + 3]);
                double z = Double.parseDouble(elems[i + 4]);
                int nrSamples = Integer.parseInt(elems[i + 5]);

                // extrapolate
//                System.out.println(origP + "\t" + z);
                double correlation3 = ZScores.zScoreToCorrelation(z, nrSamples);
                double r2 = correlation3 * correlation3;
                ln += "\t" + r2;
//                double newZ = ZScores.extrapolateZScore(nrSamples, nrNewSamples, z);
//                double newP = ZScores.zToP(newZ);

//                System.out.println("p " + origP + "\tz " + z + "\tn " + nrSamples + "\tnewz " + newZ + "\tnewp " + newP + "\tc " + correlation3);
//                row += "\t" + newP + "\t" + newZ;


            }
            toAdd.add(ln);
            probeToLn.put(probe, ctr);
            ctr++;
//            System.exit(0);
            elems = tf.readLineElems(TextFile.tab);
        }

        tf.close();

        ProbeTranslation pb = new ProbeTranslation();
        HashMap<String, String> ht12toMeta = pb.getProbeTranslation(probeTranslationFile, "HumanHT-12_V3_0_R2_11283641_A.txt", "Probe");




        TextFile tfo = new TextFile(outfile, TextFile.W);
        TextFile tfe = new TextFile(transEQTLFile, TextFile.R);
        String[] telems = tfe.readLineElems(TextFile.tab);

        while (telems.length == 1) {
            tfo.writeln(telems[0]);
            telems = tfe.readLineElems(TextFile.tab);
        }

        String header = Strings.concat(telems, Strings.tab) + "\tTransEQTLR2" + newHeader + "\tEGCUT Uncorrected R2\tEGCUT Corrected R2";
        tfo.writeln(header);
        telems = tfe.readLineElems(TextFile.tab);
        while (telems != null) {
            String snp = telems[1];
            String probe = telems[4];
            double z = Double.parseDouble(telems[9]);
            String[] zscoreStr = telems[10].split(";");
            int totalSamples = 0;
            for (int q = 0; q < zscoreStr.length; q++) {
                try {
                    Double dsz = Double.parseDouble(zscoreStr[q]);
                    totalSamples += sampleSizes[q];
                } catch (NumberFormatException e) {
                    // nothing to do here...
                }
            }

            double correlation3 = ZScores.zScoreToCorrelation(z, totalSamples);
            double r2 = correlation3 * correlation3;
            String probeMetaId = ht12toMeta.get(probe);

            Integer correspondingendos = probeToLn.get(probeMetaId);

            String zCorr = correctedZSCore.get(snp + "-" + probe);
            String zUnCorr = uncorrectedZSCore.get(snp + "-" + probe);
            if (zCorr == null) {
                zCorr = "N/A";
            }
            if (zUnCorr == null) {
                zUnCorr = "N/A";
            }
            String lnOut = Strings.concat(telems, Strings.tab) + "\t" + r2 + toAdd.get(correspondingendos) + "\t" + zUnCorr + "\t" + zCorr;
            tfo.writeln(lnOut);
            telems = tfe.readLineElems(TextFile.tab);
        }
        tfe.close();
        tfo.close();
    }

    public void test() {
        Correlation.correlationToZScore(2573);
        System.out.println(Correlation.convertCorrelationToZScore(2573, 0.160934769));
    }
}
