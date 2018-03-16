/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Gene;
import umcg.genetica.ensembl.Features;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class DetermineDistanceBetweenSNPandTSS {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String ensemblAnnotation = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/All.0.96%Identity-Merged-PerProbe-UniqueMappings-Ensembl70.txt-FullEnsemblAnnotation.txt";
        String ensemblFeatures = "/Volumes/iSnackHD/Data/GenomeSequences/Ensembl70_HG19/structures/structures_b70.txt";
        String eqtlfile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/2014-03-DistanceToTSS/EQTLs.txt";
        String snpfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12v3-1KGImputed/SNPMappings.txt";

        try {
            Features featureset = new Features();
            featureset.loadAnnotation(ensemblFeatures);

            HashMap<String, Integer> probeToTSS = new HashMap<String, Integer>();

            TextFile tf = new TextFile(ensemblAnnotation, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String probe = elems[5];

                String gene = elems[6];

                String[] geneElems = gene.split(";");
                for (String geneElem : geneElems) {
                    Gene g = featureset.getGeneHash().get(geneElem);
                    if (g != null) {
                        int start = g.getStart();
                        if (g.getStrand() < 0) {
                            start = g.getEnd();
                        }
                        probeToTSS.put(probe, start);
                    }
                }

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            HashSet<String> snps = new HashSet<String>();
            TextFile tf2 = new TextFile(eqtlfile, TextFile.R);
            tf2.readLine();
            String[] elems2 = tf2.readLineElems(TextFile.tab);
            while (elems2 != null) {
                String snp = elems2[1];
                snps.add(snp);
                elems2 = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

            TextFile tf3 = new TextFile(snpfile, TextFile.R);
            String[] elems3 = tf3.readLineElems(TextFile.tab);
            HashMap<String, Integer> snpToPos = new HashMap<String, Integer>();
            while (elems3 != null) {
                int pos = Integer.parseInt(elems3[1]);
                String snp = elems3[2];
                if (snps.contains(snp)) {
                    snpToPos.put(snp, pos);
                }
                elems3 = tf3.readLineElems(TextFile.tab);
            }
            tf3.close();

            tf2.open();
            TextFile tf4 = new TextFile(eqtlfile + "WithDistance.txt", TextFile.W);
            tf4.writeln(tf2.readLine() + "\tSNPPos\tGeneTSS");
            elems2 = tf2.readLineElems(TextFile.tab);

            while (elems2 != null) {
                String snp = elems2[1];
                String probe = elems2[4];

                Integer pos = snpToPos.get(snp);
                Integer pos2 = probeToTSS.get(probe);

                tf4.writeln(Strings.concat(elems2, Strings.tab) + "\t" + pos + "\t" + pos2);

                elems2 = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
