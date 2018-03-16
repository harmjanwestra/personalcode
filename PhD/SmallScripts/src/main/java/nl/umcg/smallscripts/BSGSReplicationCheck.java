/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class BSGSReplicationCheck {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        BSGSReplicationCheck s = new BSGSReplicationCheck();
        s.run();

    }

    public void run() {
        // TODO code application logic here
//         String bsgsfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/BSGS_trans_eQTL_replication_results2.csv";
        String bsgsfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/ParsedSorted/eQTLs-HT12v3.txt";
//        String bsgsfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/ParsedSortedFilteredForFDR0.05InMeta/eQTLsMinorAlleleInconsistenciesCorrected.txt";
        String ttgenotypeds = "/Data/GeneticalGenomicsDatasets/BloodHT12ImputeTriTyper/";
        double threshold = 0.95;

        try {
            HashSet<String> uniqueSNPs = new HashSet<String>();

            TextFile tf = new TextFile(bsgsfile, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                uniqueSNPs.add(elems[1]);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TriTyperGenotypeData ds = new TriTyperGenotypeData(ttgenotypeds);
            SNPLoader loader = ds.createSNPLoader();

            DetermineLD ldcalc = new DetermineLD();

            String[] snpsInFile = uniqueSNPs.toArray(new String[0]);

            HashSet<Pair<String, String>> linkedPairs = new HashSet<Pair<String, String>>();
            System.out.println(snpsInFile.length + " SNPs in file");
            ProgressBar pb = new ProgressBar(snpsInFile.length);


            HashMap<Pair<String, String>, Double> r2PerPair = new HashMap<Pair<String, String>, Double>();
            for (int i = 0; i < snpsInFile.length; i++) {
                Integer snpId1 = ds.getSnpToSNPId().get(snpsInFile[i]);
                if (snpId1 != null) {
                    SNP snpObj1 = ds.getSNPObject(snpId1);
                    loader.loadGenotypes(snpObj1);
                    for (int j = i + 1; j < snpsInFile.length; j++) {
                        Integer snpId2 = ds.getSnpToSNPId().get(snpsInFile[j]);

                        if (snpId2 != null && ds.getChr(snpId1) == ds.getChr(snpId2)) {
                            SNP snpObj2 = ds.getSNPObject(snpId2);
                            loader.loadGenotypes(snpObj2);

                            double r2 = ldcalc.getRSquared(snpObj1, snpObj2, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);

                            if (r2 >= threshold) {


                                Pair<String, String> p = new Pair<String, String>(snpsInFile[i], snpsInFile[j]);
                                Pair<String, String> p2 = new Pair<String, String>(snpsInFile[j], snpsInFile[i]);

                                if (!linkedPairs.contains(p) && !linkedPairs.contains(p2)) {
                                    linkedPairs.add(p);
                                    r2PerPair.put(p, r2);
                                }

//                                System.out.println(p.toString() + "\t" + r2);

                            }

                            snpObj2.clearGenotypes();
                        }

                    }
                    snpObj1.clearGenotypes();
                }
                pb.iterate();
            }
            pb.close();

            System.out.println(linkedPairs.size() + " linked pairs of SNPs");

            System.out.println("");
            int pairswithissues = 0;
            for (Pair<String, String> p : linkedPairs) {
                tf.open();
                tf.readLine();
                Double r2 = r2PerPair.get(p);
//                System.out.println(p.toString() + "\t" + r2);
                elems = tf.readLineElems(TextFile.tab);




                HashSet<String> uniqueProbes = new HashSet<String>();
                HashMap<String, Container> data = new HashMap<String, Container>();

                while (elems != null) {
                    String snp = elems[1];
                    if (snp.equals(p.getLeft()) || snp.equals(p.getRight())) {
//                        System.out.println(Strings.concat(elems, Strings.tab, 0, 8));
                        String probe = elems[4];
                        double z = Double.parseDouble(elems[10]);
                        String alleleAssessed = elems[9];
//                        String minor = elems[5];


                        Container c = data.get(probe);
                        if (c == null) {
                            c = new Container();
                        }

                        if (snp.equals(p.getLeft())) {
                            c.assesL = alleleAssessed;
                            c.zL = z;
//                            c.minorL = minor;
                            c.probe = probe;
                            c.snpL = snp;
                        } else {
                            c.assesR = alleleAssessed;
                            c.zR = z;
//                            c.minorR = minor;
                            c.probe = probe;
                            c.snpR = snp;
                        }
                        uniqueProbes.add(probe);
                        data.put(probe, c);

                    }
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();

                String[] allProbes = uniqueProbes.toArray(new String[0]);

                Container firstFx = null;
                String alleleL = null;
                String alleleR = null;
                String minorL = null;
                String minorR = null;

                boolean isFirstEffect = false;
                boolean firstEffectHasOppositeDirection = false;

                boolean hasIssues = false;
                for (String probe : allProbes) {
                    // irrespective of the direction, set initial alleles first...
                    Container c = data.get(probe);

                    if (c.assesL != null && c.assesR != null) {
                        if (alleleL == null) {
                            alleleL = c.assesL;
                            alleleR = c.assesR;

                            isFirstEffect = true;
                            if ((c.zL >= 0 && c.zR >= 0) || (c.zL < 0 && c.zR < 0)) {
                                firstEffectHasOppositeDirection = false;
                            } else {
                                firstEffectHasOppositeDirection = true;
                            }
                            firstFx = c;
//                            System.out.println(c.toString());
                        } else {
                            isFirstEffect = false;
                            // first check whether the assessed alleles have stayed identical...

                            if (alleleL.equals(c.assesL) && alleleR.equals(c.assesR)) {
                                // the alleles have stayed the same.. Is the direction of effects also identical?
//                                System.out.println(c.toString());
                                if ((c.zL >= 0 && c.zR >= 0) || (c.zL < 0 && c.zR < 0)) { // effect has identical direction between SNPs
                                    if (firstEffectHasOppositeDirection) { // alleles haven't changed, but they used to show an opposite direction in the first effect
                                        System.out.println("First:\t" + firstFx.toString() + "\t" + r2);
                                        System.out.println("Second:\t" + c.toString() + "\t" + r2);
//                                        System.out.println("");
                                        hasIssues = true;
                                    }
                                } else { // effect has opposite direction between SNPs
                                    if (!firstEffectHasOppositeDirection) { // the first effect did not have an opposite direction
                                        System.out.println("First:\t" + firstFx.toString() + "\t" + r2);
                                        System.out.println("Second:\t" + c.toString() + "\t" + r2);
//                                        System.out.println("");
                                        hasIssues = true;
                                    }
                                }
                            } else { // the assessed alleles have changed
                                if (!alleleL.equals(c.assesL) && !alleleR.equals(c.assesR)) {
                                } else {
                                    if ((c.zL >= 0 && c.zR >= 0) || (c.zL < 0 && c.zR < 0)) {  // effect has identical direction between SNPs
                                        if (!firstEffectHasOppositeDirection) { // the direction is now identical, although the alleles have changed. Was the direcgtion of the effect also identical with the other alleles?
                                            System.out.println("First:\t" + firstFx.toString() + "\t" + r2);
                                            System.out.println("Second:\t" + c.toString() + "\t" + r2);
//                                            System.out.println("");
                                            hasIssues = true;
                                        }
                                    } else {
                                        if (firstEffectHasOppositeDirection) { // the direction is now opposite, although the alleles have changed. Was the direcgtion of the effect also opposite with the other alleles?
                                            System.out.println("First:\t" + firstFx.toString() + "\t" + r2);
                                            System.out.println("Second:\t" + c.toString() + "\t" + r2);
//                                            System.out.println("");
                                            hasIssues = true;
                                        }
                                    }
                                }
                            }
                        }
                    }


                }
//                System.out.println("");
                if (hasIssues) {
                    pairswithissues++;
                }


            }

            System.out.println(pairswithissues + " SNP pairs have issues.");

        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public class Container {
        // start stepping through the array from the beginning

        public String minorR = null;
        public String minorL = null;
        public String assesR = null;
        public String assesL = null;
        double zR = 0;
        double zL = 0;
        String snpR = null;
        String snpL = null;
        String probe = null;

        @Override
        public String toString() {
            return probe + "\t" + snpR + "\t" + zR + "\t" + assesR + "\t"
                    + snpL + "\t" + zL + "\t" + assesL;
        }
    }
}
