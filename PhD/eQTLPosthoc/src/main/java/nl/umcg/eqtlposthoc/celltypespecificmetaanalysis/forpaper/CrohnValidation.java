/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CrohnValidation {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String celltypevector = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysis/Vector-CellTypeInteractionZScore.txt";
            String binomialTestFile = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-12-04-2-NoGWASPValueFilter-BinomialDiseaseTestFDR0.05/FisherExactOutput-[CellTypeInteractionZScore]-Pruning0.2-BgBuildWithGWASSNPsOnlyfalse.txt";
            String replicationFile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/2013-11-ManuscriptNatureMethods-Draft1/Supplementary/Supplementary Table 2 - Replication Effect sizes and mean expression.txt";
            String output = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-12-04-2-NoGWASPValueFilter-BinomialDiseaseTestFDR0.05/ListOfEffectsForCrohnsDisease.txt";

            String querytrait = "Crohn's disease";

            TextFile tf2 = new TextFile(binomialTestFile, TextFile.R);
            String[] elems = tf2.readLineElems(TextFile.tab);

            String[] prunedSNPs = null;
            
            while (elems != null) {
                String trait = elems[5];
                if (trait.equals(querytrait)) {
                    prunedSNPs = elems[elems.length - 1].split(",");
                }
                elems = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

            if (prunedSNPs == null) {
                System.err.println("No SNPs found");
                System.exit(-1);
            }

            HashSet<String> hashTraitSNPs = new HashSet<String>();
            hashTraitSNPs.addAll(Arrays.asList(prunedSNPs));
            TextFile tf = new TextFile(celltypevector, TextFile.R);
            tf.readLine();
            elems = tf.readLineElems(TextFile.tab);
            HashMap<String, Pair<String, Double>> traitEQTLs = new HashMap<String, Pair<String, Double>>();
            HashSet<String> selectedEQTLs = new HashSet<String>();
            while (elems != null) {
                if (elems.length > 1) {
                    String eqtl = elems[0];
                    String dbl = elems[1];
                    Double dbl2 = Double.parseDouble(dbl);

                    String snp = eqtl.split("-")[0];
                    if (hashTraitSNPs.contains(snp)) {
                        traitEQTLs.put(snp, new Pair<String, Double>(eqtl, dbl2));
                        selectedEQTLs.add(eqtl);
                    }
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

//            ZScores.
            
            TextFile tf3 = new TextFile(replicationFile, TextFile.R);
            TextFile tf4 = new TextFile(output, TextFile.W);

            tf4.writeln(tf3.readLine());
            elems = tf3.readLineElems(TextFile.tab);
            while (elems != null) {
                String eqtl = elems[0] + "-" + elems[1];
                if (selectedEQTLs.contains(eqtl)) {
                    tf4.writeln(Strings.concat(elems, Strings.tab));
                }

                elems = tf3.readLineElems(TextFile.tab);
            }

            tf4.close();
            tf3.close();

        } catch (IOException e) {
            e.printStackTrace();

        }
    }

}
