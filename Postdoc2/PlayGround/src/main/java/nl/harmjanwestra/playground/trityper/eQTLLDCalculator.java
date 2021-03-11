package nl.harmjanwestra.playground.trityper;

import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class eQTLLDCalculator {

    public static void main(String[] args) {
        eQTLLDCalculator e = new eQTLLDCalculator();
        System.out.println("eQTL file top association LD calculator.");
        if (args.length < 6) {

            System.out.println("efile1 efile 2 trityper individuals.txt probeannotation.txt outputfile");
        } else {
            try {
                e.run(args[0], args[1], args[2], args[3], args[4], args[5]);
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        }
    }

    public void run(String efile1, String efile2, String triTyperDataset, String individualSubset, String geneAnnotation, String outputFile) throws IOException {

        HashMap<String, String> geneChrs = new HashMap<>();
        HashMap<String, Integer> genePos = new HashMap<>();
        System.out.println("Loading gene annotation: " + geneAnnotation);
        TextFile tfa = new TextFile(geneAnnotation, TextFile.R);
        tfa.readLine();
        String[] elems = tfa.readLineElems(TextFile.tab);
        while (elems != null) {
            String gene = elems[1];
            String chr = elems[3];
            int pos = Integer.parseInt(elems[4]);
            geneChrs.put(gene, chr);
            genePos.put(gene, pos);
            elems = tfa.readLineElems(TextFile.tab);
        }
        tfa.close();

        Triple<HashMap<String, String>, HashMap<String, String>, HashMap<String, String>> topsnpq1 = readTopSNPs(efile1);
        HashMap<String, String> topsnp1 = topsnpq1.getLeft();
        HashMap<String, String> topsnp1chr = topsnpq1.getMiddle();
        HashMap<String, String> topsnp1pos = topsnpq1.getRight();

        Triple<HashMap<String, String>, HashMap<String, String>, HashMap<String, String>> topsnpq2 = readTopSNPs(efile2);
        HashMap<String, String> topsnp2 = topsnpq2.getLeft();
        HashMap<String, String> topsnp2chr = topsnpq2.getMiddle();
        HashMap<String, String> topsnp2pos = topsnpq2.getRight();

        // preload genotypes
        System.out.println("Loading genotypes");
        LDCalculator ldcalc = new LDCalculator();
        TriTyperGenotypeData data = new TriTyperGenotypeData(triTyperDataset);
        if (individualSubset != null) {
            HashSet<String> indsToInclude = ldcalc.loadSet(individualSubset);
            System.out.println("Filtering individuals using: " + individualSubset + " having " + indsToInclude.size() + " individuals");
            String[] individuals = data.getIndividuals();
            Boolean[] include = new Boolean[individuals.length];
            for (int i = 0; i < individuals.length; i++) {
                if (indsToInclude.contains(individuals[i])) {
                    include[i] = true;
                } else {
                    include[i] = false;
                }
            }
            data.setIsIncluded(include);
        }


        ArrayList<String> allSNPs = new ArrayList<String>();
        allSNPs.addAll(topsnp1.values());
        allSNPs.addAll(topsnp2.values());
        System.out.println(allSNPs.size() + " total snps to load");
        HashMap<String, SNP> genotypes = ldcalc.loadSNPs(allSNPs, data, false);
        System.out.println(genotypes.size() + " SNP genotypes read");

        TextFile output = new TextFile(outputFile, TextFile.W);
        output.writeln("Gene\tGeneChr\tGenePos\tSNP1\tPos1\tTSSDistSNP1\tSNP2\tPos2\tTSSDistSNP2\tSNPSNPDistance\tRsquared\tDprime");
        for (String gene : topsnp1.keySet()) {
            String snp1 = topsnp1.get(gene);
            String snp2 = topsnp2.get(gene);

            String genechr = geneChrs.get(gene);
            int genepos = genePos.get(gene);

            if (snp1 != null && snp2 != null) {

                int pos1 = Integer.parseInt(topsnp1pos.get(gene));
                int pos2 = Integer.parseInt(topsnp2pos.get(gene));

                int snpdist = Math.abs(pos1 - pos2);
                int tssdist1 = pos1 - genepos;
                int tssdist2 = pos2 - genepos;
                Double rsquared = null;
                Double dprime = null;

                if (snp1.equals(snp2)) {
                    rsquared = 1d;
                    dprime = 1d;
                } else {
                    SNP obj1 = genotypes.get(snp1);
                    SNP obj2 = genotypes.get(snp2);
                    if (obj1 != null && obj2 != null) {
                        DetermineLD ld = new DetermineLD();
                        Pair<Double, Double> ldvals = ld.getLD(obj1, obj2, data, 1, false);
                        rsquared = ldvals.getRight();
                        dprime = ldvals.getLeft();
                    }
                }

                String rsquaredStr = "" + rsquared;
                String dprimeStr = "" + dprime;
                if (rsquared == null) {
                    rsquaredStr = "-";
                    dprimeStr = "-";
                }

                output.writeln(gene
                        + "\t" + genechr
                        + "\t" + genepos
                        + "\t" + snp1
                        + "\t" + pos1
                        + "\t" + tssdist1
                        + "\t" + snp2
                        + "\t" + pos2
                        + "\t" + tssdist2
                        + "\t" + snpdist + "\t" + rsquaredStr + "\t" + dprimeStr);
            }
        }
        output.close();

    }

    private Triple<HashMap<String, String>, HashMap<String, String>, HashMap<String, String>> readTopSNPs(String efile1) throws IOException {
        HashMap<String, String> topsnp = new HashMap<>();
        System.out.println("Reading: " + efile1);
        TextFile tf = new TextFile(efile1, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        HashMap<String, Double> lowestP = new HashMap<>();
        HashMap<String, String> lowestPChr = new HashMap<>();
        HashMap<String, String> lowestPPos = new HashMap<>();
        while (elems != null) {
            double p = Double.parseDouble(elems[0]);
            String snp = elems[1];
            String chr = elems[2];
            String pos = elems[3];
            String gene = elems[4];
            Double q = lowestP.get(gene);
            if (q == null || p < q) {
                topsnp.put(gene, snp);
                lowestP.put(gene, p);
                lowestPChr.put(gene, chr);
                lowestPPos.put(gene, pos);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(topsnp.size() + " top eQTLs");
        return new Triple<HashMap<String, String>, HashMap<String, String>, HashMap<String, String>>(topsnp, lowestPChr, lowestPPos);
    }

}
