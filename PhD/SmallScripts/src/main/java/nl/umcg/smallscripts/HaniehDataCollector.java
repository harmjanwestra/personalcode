/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class HaniehDataCollector {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            HaniehDataCollector d = new HaniehDataCollector();
            
            String eqtlfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLs.txt-FilteredForProbeLevelFDR.txt.gz";
            String outfile = "/Volumes/iSnackHD/Data/Projects/2013-IVAnalysis/2013-01-08-SNPsWithSingleCisEQTL.txt";
            double threshold = 1;
            d.run(eqtlfile, threshold, outfile);
            
            String outfile2 = "/Volumes/iSnackHD/Data/Projects/2013-IVAnalysis/2013-01-08-SNPsWithSingleCisEQTL-FilteredForTopEffects.txt";
            d.filterTopEffects(eqtlfile, outfile, outfile2);
            
            String outfile3 = "/Volumes/iSnackHD/Data/Projects/2013-IVAnalysis/2013-01-08-SNPsWithSingleCisEQTL-FilteredForTopEffects-WithAverageGeneExpression.txt";
            String outfile4 = "/Volumes/iSnackHD/Data/Projects/2013-IVAnalysis/2013-01-08-AverageGeneExpressionOfOtherProbes.txt";
            String expressionFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/eQTLsVsGeneExpression/cisVsGeneExpression/cisProbeLevelFDR0.05-MeanExpressionOfProbes.txt";
            d.combineWithAverageGeneExpression(outfile2, expressionFile, outfile3, outfile4);
            
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String eqtlfile, double threshold, String outfile) throws IOException {
        // read eQTL file.. process until threshold has been met

        TextFile tf = new TextFile(eqtlfile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        HashMap<String, HashSet<String>> snpToGeneMap = new HashMap<String, HashSet<String>>();
        HashSet<String> snps = new HashSet<String>();
        HashMap<String, String> probeToGeneMap = new HashMap<String, String>();
        while (elems != null) {

            double p = Double.parseDouble(elems[0]);
            if (p > threshold) {

                break;
            } else {
                String snp = elems[1];
                String probe = elems[4];
                String gene = elems[eQTLTextFile.HUGO];

                // select snps with a single gene eQTL
                HashSet<String> probes = snpToGeneMap.get(snp);
                if (probes == null) {
                    probes = new HashSet<String>();
                }
                probes.add(probe);
                snpToGeneMap.put(snp, probes);
                snps.add(snp);
                probeToGeneMap.put(probe, gene);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile out = new TextFile(outfile, TextFile.W);
        for (String snp : snps) {
            HashSet<String> probesForSNP = snpToGeneMap.get(snp);
            if (probesForSNP.size() == 1) {
                for (String probe : probesForSNP) {
                    out.writeln(snp + "\t" + probe + "\t" + probeToGeneMap.get(probe));
                }
            } else {

                HashSet<String> genes = new HashSet<String>();
                for (String probe : probesForSNP) {
                    String gene = probeToGeneMap.get(probe);
                    if (gene.trim().length() == 0) {
                        // no name for gene
                        genes.add("-");
                    } else {
                        genes.add(gene);
                    }
                }

                if (genes.size() == 1) {
                    for (String probe : probesForSNP) {
                        String gene = probeToGeneMap.get(probe);
                        if (gene.trim().length() == 0) {
                            out.writeln(snp + "\t" + probe + "\t-");
                        } else {
                            out.writeln(snp + "\t" + probe + "\t" + probeToGeneMap.get(probe));
                        }
                    }
                }

            }

        }
        out.close();
    }

    private void filterTopEffects(String eqtlfile, String outfile, String outfile2) throws IOException {

        HashSet<Pair<String, String>> snpProbePairsAllowed = new HashSet<Pair<String, String>>();
        TextFile in1 = new TextFile(outfile, TextFile.R);
        String[] elems = in1.readLineElems(TextFile.tab);
        while (elems != null) {
            snpProbePairsAllowed.add(new Pair<String, String>(elems[0], elems[1]));
            elems = in1.readLineElems(TextFile.tab);
        }
        in1.close();

        HashSet<String> probesVisited = new HashSet<String>();

        TextFile tf = new TextFile(eqtlfile, TextFile.R);
        TextFile tfout = new TextFile(outfile2, TextFile.W);
        String header = tf.readLine();
        tfout.writeln(header);
        String[] eelems = tf.readLineElems(TextFile.tab);
        while (eelems != null) {

            String snp = eelems[1];
            String probe = eelems[4];
            if (snpProbePairsAllowed.contains(new Pair<String, String>(snp, probe)) && !probesVisited.contains(probe)) {
                tfout.writeln(Strings.concat(eelems, Strings.tab));
                probesVisited.add(probe);
            }
            eelems = tf.readLineElems(TextFile.tab);
        }
        tfout.close();
        tf.close();
    }

    private void combineWithAverageGeneExpression(String outfile2, String expressionFile, String outfile3, String outfile4) throws IOException {
        HashMap<String, String> probeToAnnotation = new HashMap<String, String>();
        TextFile tf = new TextFile(expressionFile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);

        // nr	metaId	ht12Id	genename	isEQTL	meanExp	stdev	Z	CV	metaZ
        // 0	17272	6450255	MACC1	    false	6.562089831982889	0.08387683290521249	-0.4788440298858032	0.012782030580624822	null

        while (elems != null) {
            if (elems.length > 9) {
                String probe = elems[1];
                String out = elems[2] + "\t" + elems[5] + "\t" + elems[6] + "\t" + elems[7] + "\t" + elems[8];
                probeToAnnotation.put(probe, out);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tfout = new TextFile(outfile3, TextFile.W);
        TextFile tfin = new TextFile(outfile2, TextFile.R);
        elems = tfin.readLineElems(TextFile.tab); // skip header
        elems = tfin.readLineElems(TextFile.tab);

        HashSet<String> selectedProbes = new HashSet<String>();
        tfout.writeln("PValue\tSNP\tProbe\tGene\tEQTLMetaZ\tHT12v3id\tMeanExp\tStDev\tExpZ\tExpCV");
        while (elems != null) {

            selectedProbes.add(elems[4]);
            tfout.writeln(elems[0] + "\t" + elems[1] + "\t" + elems[4] + "\t" + elems[eQTLTextFile.HUGO] + "\t" + elems[eQTLTextFile.METAZ] + "\t" + probeToAnnotation.get(elems[4]));

            elems = tfin.readLineElems(TextFile.tab);
        }
        tfin.close();
        tfout.close();

        TextFile out2 = new TextFile(outfile4, TextFile.W);
        tf.open();
        tf.readLineElems(TextFile.tab);
        tf.readLineElems(TextFile.tab);
        elems = tf.readLineElems(TextFile.tab);
        out2.writeln("MetaID\tHT12v3id\tMeanExp\tStDev\tExpZ\tExpCV");
        while (elems != null) {
            if (elems.length > 9) {
                String probe = elems[1];
                out2.writeln(elems[1] + "\t" + elems[2] + "\t" + elems[5] + "\t" + elems[6] + "\t" + elems[7] + "\t" + elems[8]);
//                probeToAnnotation.put(probe, out);
            }
            elems = tf.readLineElems(TextFile.tab);
        }

        out2.close();
        tf.close();
    }
}
