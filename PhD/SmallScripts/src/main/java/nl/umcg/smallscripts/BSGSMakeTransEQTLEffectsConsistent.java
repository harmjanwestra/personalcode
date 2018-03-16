/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class BSGSMakeTransEQTLEffectsConsistent {

    public BSGSMakeTransEQTLEffectsConsistent() {

        HashMap hashProbes = new HashMap();
        Vector vecProbes = new Vector();
        try {
            double threshold = 0.95;
            String ttgenotypeds = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/Hap2ImputedGenotypes/";
            TriTyperGenotypeData ds = new TriTyperGenotypeData(ttgenotypeds);
            String outdir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/ParsedSortedFilteredForFDR0.05InMetaMinorAlleleInconsistenciesCorrected/";
            Gpio.createDir(outdir);
            SNPLoader loader = ds.createSNPLoader();
            for (int i = 0; i < 11; i++) {

                String bsgsfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/ParsedSortedFilteredForFDR0.05InMeta/eQTLs.txt";
                String bsgsfileOut = outdir + "eQTLs.txt";
                if (i > 0) {
                    bsgsfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/ParsedSortedFilteredForFDR0.05InMeta/PermutedEQTLsPermutationRound" + i + ".txt.gz";
                    bsgsfileOut = outdir + "PermutedEQTLsPermutationRound" + i + ".txt.gz";
                }

                System.out.println(bsgsfile);
//        String ttgenotypeds = "/Data/GeneticalGenomicsDatasets/BloodHT12ImputeTriTyper/";




                DetermineLD ldcalc = new DetermineLD();


                HashSet<String> uniqueSNPs = new HashSet<String>();

                TextFile tf = new TextFile(bsgsfile, TextFile.R);
                TextFile tfOut = new TextFile(bsgsfileOut, TextFile.W);
                tfOut.writeln(tf.readLine());
                String[] elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {

                    if (!hashProbes.containsKey(elems[4])) {
                        hashProbes.put(elems[4], elems[1] + "\t" + elems[9] + "\t" + elems[10]);
                        vecProbes.add(elems[4]);
                    } else {
                        hashProbes.put(elems[4], (String) hashProbes.get(elems[4]) + "," + elems[1] + "\t" + elems[9] + "\t" + elems[10]);
                    }

                    String snp = elems[1];
                    Integer snpId = ds.getSnpToSNPId().get(snp);
                    SNP snpObj = ds.getSNPObject(snpId);
                    loader.loadGenotypes(snpObj);
                    String allele1 = BaseAnnot.toString(snpObj.getAlleles()[0]);
                    String allele2 = BaseAnnot.toString(snpObj.getAlleles()[1]);
                    String minor = BaseAnnot.toString(snpObj.getMinorAllele());

                    elems[13] = ""+862;
                    if (!elems[9].equals(minor) && snpObj.getMAF() < 0.40) {
                        System.out.println(Strings.concat(elems, Strings.tab) + "\t" + allele1 + "/" + allele2 + "\t" + minor + "\t" + elems[9].equals(minor) + "\t" + snpObj.getMAF());
                        elems[10] = String.valueOf(-Double.parseDouble(elems[10]));
                    }
                    tfOut.writeln(Strings.concat(elems, Strings.tab));

                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
                tfOut.close();

//                if (1 == 1) {
//                    System.exit(0);
//                }
//
//                System.out.println(vecProbes.size());
//
//
//                for (int p = 0; p < vecProbes.size(); p++) {
//                    String probe = (String) vecProbes.get(p);
//                    String[] eQTLs = ((String) hashProbes.get(probe)).split(",");
//                    for (int a = 0; a < eQTLs.length; a++) {
//                        String snpA = eQTLs[a].split("\t")[0];
//                        Integer snpIdA = ds.getSnpToSNPId().get(snpA);
//                        SNP snpObjA = ds.getSNPObject(snpIdA);
//                        loader.loadGenotypes(snpObjA);
//                        String allele1A = BaseAnnot.toString(snpObjA.getAlleles()[0]);
//                        String allele2A = BaseAnnot.toString(snpObjA.getAlleles()[1]);
//                        String minorA = BaseAnnot.toString(snpObjA.getMinorAllele());
//
//                        for (int b = a + 1; b < eQTLs.length; b++) {
//                            String snpB = eQTLs[b].split("\t")[0];
//                            Integer snpIdB = ds.getSnpToSNPId().get(snpB);
//                            SNP snpObjB = ds.getSNPObject(snpIdB);
//                            loader.loadGenotypes(snpObjB);
//                            String allele1B = BaseAnnot.toString(snpObjB.getAlleles()[0]);
//                            String allele2B = BaseAnnot.toString(snpObjB.getAlleles()[1]);
//                            String minorB = BaseAnnot.toString(snpObjB.getMinorAllele());
//
//                            double r2 = ldcalc.getRSquared(snpObjA, snpObjB, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
//
//                            if (eQTLs[a].split("\t")[1].equals(minorA) && eQTLs[b].split("\t")[1].equals(minorB)) {
//                            } else {
//                                if (r2 > 0.8) {
//
//                                    System.out.println(probe + "\t" + eQTLs.length + "\t" + a + "\t" + b + "\t" + r2 + "\t" + eQTLs[a] + "\t" + allele1A + "/" + allele2A + "\t" + minorA + "\t" + eQTLs[b] + "\t" + allele1B + "/" + allele2B + "\t" + minorB);
//                                }
//                            }
//
//                        }
//                    }
//                }


            }
            loader.close();



        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.exit(0);
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        new BSGSMakeTransEQTLEffectsConsistent();
    }
}
