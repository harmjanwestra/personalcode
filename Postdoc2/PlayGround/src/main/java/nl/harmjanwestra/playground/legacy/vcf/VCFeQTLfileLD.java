package nl.harmjanwestra.playground.legacy.vcf;

import nl.harmjanwestra.playground.trityper.LDCalculator;
import nl.harmjanwestra.playground.trityper.eQTLLDCalculator;
import org.checkerframework.checker.units.qual.A;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.features.SNPFeature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.ChrAnnotation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class VCFeQTLfileLD {


    public static void main(String[] args) {
        VCFeQTLfileLD e = new VCFeQTLfileLD();
        System.out.println("eQTL file top association LD calculator.");
        if (args.length < 6) {

            System.out.println("efile1 efile2 vcf individuals.txt probeannotation.txt outputfile");
        } else {
            try {
                e.run(args[0], args[1], args[2], args[3], args[4], args[5]);
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        }
    }

    public void run(String efile1, String efile2, String vcftemplate, String individualSubset, String geneAnnotation, String outputFile) throws IOException {

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


        ArrayList<String> allSNPs = new ArrayList<String>();
        allSNPs.addAll(topsnp1.values());
        allSNPs.addAll(topsnp2.values());
        System.out.println(allSNPs.size() + " total snps to load");

        // sort by chr
        ArrayList<Pair<String, Integer>> genep = new ArrayList<>();
        for (String gene : topsnp1.keySet()) {
            Integer chr = Chromosome.parseChr(geneChrs.get(gene)).getNumber();
            genep.add(new Pair<>(gene, chr, Pair.SORTBY.RIGHT));
        }

        Collections.sort(genep);

        TextFile output = new TextFile(outputFile, TextFile.W);
        output.writeln("Gene\tGeneChr\tGenePos\tSNP1\tPos1\tTSSDistSNP1\tSNP2\tPos2\tTSSDistSNP2\tSNPSNPDistance\tRsquared\tDprime");
        ProgressBar pb = new ProgressBar(topsnp1.keySet().size(), "Calculating LD");
        Integer prevChr = null;
        VCFTabix tabix = null;
        boolean[] samplefilter = null;

        for (Pair<String, Integer> genepair : genep) {
            String gene = genepair.getLeft();
            Integer chr = genepair.getRight();

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
                    int geneChr = ChrAnnotation.parseChr(genechr);
                    if (prevChr == null || prevChr != geneChr) {
                        if (tabix != null) {
                            tabix.close();
                        }

                        String vcf = vcftemplate.replaceAll("CHR", "" + geneChr);
                        System.out.println("Opening: " + vcf);
                        tabix = new VCFTabix(vcf);
                        samplefilter = tabix.getSampleFilter(individualSubset);
                        prevChr = geneChr;
                    }

                    VCFVariant var1 = tabix.getVariant(asFeature(geneChr, pos1), samplefilter);
                    VCFVariant var2 = tabix.getVariant(asFeature(geneChr, pos2), samplefilter);
                    if (var1 != null && var2 != null) {
                        DetermineLD ld = new DetermineLD();
                        Pair<Double, Double> ldvals = ld.getLD(var1, var2);
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
            pb.iterate();
        }
        if (tabix != null) {
            tabix.close();
        }
        pb.close();
        output.close();

    }

    private Feature asFeature(int geneChr, int pos1) {
        SNPFeature feature = new SNPFeature(Chromosome.parseChr("" + geneChr), pos1, pos1);
        return feature;
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
