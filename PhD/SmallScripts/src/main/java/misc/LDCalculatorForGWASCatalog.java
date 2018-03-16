/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class LDCalculatorForGWASCatalog {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String snpsTested = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysis/Vector-CellTypeInteractionZScore.txt";
        String gwasCatalog = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-09-23-gwascatalog-WithAllIBDLoci.txt";
        String gwasTrait = "Nature2012-IBD+Crohn-NoUCSpecific";
        double significancethreshold = 1;
        double ldthreshold = 0.98;
        String thousandgenomes = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap2r24-CEU/";//"/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/Merged/";//  //  //;

        try {
            GWASCatalog catalog = new GWASCatalog(gwasCatalog, significancethreshold);

            GWASSNP[] snps = catalog.getSNPsForTraitContainingKey(gwasTrait);
            String[] traitSNPs = new String[snps.length];
            for (int i = 0; i < traitSNPs.length; i++) {
                traitSNPs[i] = snps[i].getName();
            }

            HashSet<String> uniqueSNPsInteractionTerm = new HashSet<String>();
            TextFile tf = new TextFile(snpsTested, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                uniqueSNPsInteractionTerm.add(elems[0].split("-")[0]);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            System.out.println(uniqueSNPsInteractionTerm.size() + " SNPs loaded from vector");

            TriTyperGenotypeData ds = new TriTyperGenotypeData(thousandgenomes);
            System.out.println("Done loading data..");
            SNPLoader loader = ds.createSNPLoader();
            DetermineLD ld = new DetermineLD();

            for (int i = 0; i < traitSNPs.length; i++) {
                String traitSNP = traitSNPs[i];
                if (!uniqueSNPsInteractionTerm.contains(traitSNP)) {
                    Integer snpId = ds.getSnpToSNPId().get(traitSNP);
                    if (snpId != null) {
                        SNP traitSNPObj = ds.getSNPObject(snpId);
                        loader.loadGenotypes(traitSNPObj);
                        boolean hasProxy = false;
                        for (String interactionSNP : uniqueSNPsInteractionTerm) {
                            Integer snpId2 = ds.getSnpToSNPId().get(interactionSNP);
                            if (snpId2 != null) {
                                SNP interactionSNPObj = ds.getSNPObject(snpId2);
                                if (interactionSNPObj.getChr() == traitSNPObj.getChr() && traitSNPObj.getChr() > 0) {
                                    loader.loadGenotypes(interactionSNPObj);

                                    Pair<Double, Double> pair = ld.getLD(traitSNPObj, interactionSNPObj, ds, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                                    double r2 = pair.getRight();
                                    double dp = pair.getLeft();
                                    if (r2 >= ldthreshold) {
                                        System.out.println(traitSNP + "\t" + interactionSNP + "\t" + r2 + "\t" + traitSNPObj.getChr() + "\t" + interactionSNPObj.getChr() + "\t" + traitSNPObj.getChrPos() + "\t" + interactionSNPObj.getChrPos() + "\t" + Math.abs(traitSNPObj.getChrPos() - interactionSNPObj.getChrPos()));
                                        hasProxy = true;
                                    }
                                    interactionSNPObj.clearGenotypes();
                                }

                            }
                        }

                        if (!hasProxy) {
                            System.out.println(traitSNP + "\tnull\thas no proxy");
                        }

                        traitSNPObj.clearGenotypes();
                    } else {
                        System.out.println(traitSNP + "\tnull\tnot present in ref");
                    }
                } else {
                    System.out.println(traitSNP + "\t" + traitSNP + "\tSNP already in list");
                }

            }

            loader.close();

        } catch (IOException e) {
            e.printStackTrace();

        }

    }

}
