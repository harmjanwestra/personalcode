package nl.harmjanwestra.playground.biogen.genotypes;

import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import nl.harmjanwestra.playground.legacy.vcf.filter.variantfilters.VCFVariantCallRateFilter;
import nl.harmjanwestra.playground.legacy.vcf.filter.variantfilters.VCFVariantFilters;
import nl.harmjanwestra.playground.legacy.vcf.filter.variantfilters.VCFVariantHWEPFilter;
import nl.harmjanwestra.playground.legacy.vcf.filter.variantfilters.VCFVariantMAFFilter;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.IntStream;

public class CalculateGRM {

    public static void main(String[] args) {

        if (args.length < 2) {
            System.out.println("Usage: inputCHR.vcf.gz output [variantlimit.txt]");
        } else {
            CalculateGRM s = new CalculateGRM();
            try {
                String varlimit = null;
                if (args.length > 2) {
                    varlimit = args[2];
                }
                s.run(args[0], args[1], varlimit);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public void run(String input, String output, String varlimit) throws IOException {

        HashSet<String> allowedVariants = null;
        if (varlimit != null) {
            allowedVariants = new HashSet<>();

            TextFile tf = new TextFile(varlimit, TextFile.R);
            String ln = tf.readLine();
            while (ln != null) {
                allowedVariants.add(ln);
                ln = tf.readLine();
            }
            tf.close();
            System.out.println(allowedVariants.size() + " variant IDS read from " + varlimit);
        }
        HashMap<String, Integer> sampleHash = new HashMap<>();
        ArrayList<String> samples = new ArrayList<>();

        double[][] grm = null;
        double[][] shared = null;

        int[][] total = null;

        VCFVariantFilters filters = new VCFVariantFilters();
        filters.add(new VCFVariantMAFFilter(0.10));
        filters.add(new VCFVariantCallRateFilter(0.95));
        filters.add(new VCFVariantHWEPFilter(0.0001));

        int ctr = 0;
        int processed = 0;
        int startchr = 1;
        int stopchr = 2;
        if (input.contains("CHR")) {
            startchr = 1;
            stopchr = 23;
        }
        for (int chr = startchr; chr < stopchr; chr++) {
            String vcffile = input.replace("CHR", "" + chr);
            TextFile tf = new TextFile(vcffile, TextFile.R);
            String ln = tf.readLine();
            int[] sampleOrder = null;
            while (ln != null) {
                if (ln.startsWith("#")) {
                    // header line
                    if (ln.startsWith("#CHROM")) {
                        // initialize if not yet done so
                        if (total == null) {
                            String[] elems = ln.split("\t");
                            int sctr = 0;
                            for (int i = 9; i < elems.length; i++) {
                                sampleHash.put(elems[i], sctr);
                                samples.add(elems[i]);
                                sctr++;
                            }
                            System.out.println(vcffile + " has " + samples.size() + " samples");
                            sampleOrder = new int[samples.size()];
                            for (int i = 0; i < samples.size(); i++) {
                                sampleOrder[i] = i;
                            }
                            total = new int[samples.size()][samples.size()];
                            grm = new double[samples.size()][samples.size()];
                            shared = new double[samples.size()][samples.size()];

                        } else {
                            // reorder samples to first file
                            System.out.println();
                            String[] elems = ln.split("\t");
                            sampleOrder = new int[samples.size()];
                            Arrays.fill(sampleOrder, -1);
                            System.out.println(vcffile + " has " + (elems.length - 9) + " samples, " + sampleOrder.length);
                            for (int i = 9; i < elems.length; i++) {
                                Integer id = sampleHash.get(elems[i]);
//                                System.out.println(elems[i] + " --> " + id);
//                                if (id >= samples.size()) {
//                                    System.out.println("ERRRORRRRR");
//                                    System.out.println(elems[i] + " --> " + id);
//                                    System.out.println(elems.length);
//                                    System.exit(0);
//                                }
                                if (id != null) {
                                    sampleOrder[i-9] = id;
                                }
                            }
                        }
                        System.out.println("Number of samples: " + samples.size() + "\tUnique: " + sampleHash.size());
                    }
                } else {
                    // variant line.
                    // parse variant

                    boolean parse = false;
                    if (allowedVariants != null) {
                        VCFVariant tmpvariant = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
                        if (allowedVariants.contains(tmpvariant.getId())) {
                            parse = true;
                        }
                    } else {
                        parse = true;
                    }

                    if (parse) {
                        VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
                        if (filters.passesFilters(variant)) {
                            double[] freq = variant.getAlleleFrequencies();
                            double refFreq = freq[0]; // [0.1 - 0.9]
//						System.out.println();
//						System.out.println(refFreq);

                            byte[] dosages = variant.getGenotypesAsByteVector();
                            // reorder dosages
                            byte[] tmpdosages = new byte[samples.size()];
                            byte b = -1;
                            Arrays.fill(tmpdosages, b);
                            for (int i = 0; i < sampleOrder.length; i++) {
                                if (sampleOrder[i] != -1) {
                                    tmpdosages[i] = dosages[sampleOrder[i]];
                                }
                            }

                            double[][] finalGrm = grm;
                            byte[] finaldosages = tmpdosages;
                            int[][] finalTotal = total;
                            double[][] finalShared = shared;
                            IntStream.range(0, finaldosages.length).forEach(i -> {
                                for (int j = i; j < finaldosages.length; j++) {
                                    // https://www.dissect.ed.ac.uk/documentation-grm-computation/
                                    if (finaldosages[i] != -1 && finaldosages[j] != -1) { // skip missing genotypes
                                        double a = (
                                                (Math.abs(2 - finaldosages[i]) - (2 * refFreq)) * (Math.abs(2 - finaldosages[j]) - (2 * refFreq))
                                        )
                                                / (2 * refFreq * (1 - refFreq));
                                        finalGrm[i][j] += a;
                                        if (finaldosages[i] == finaldosages[j]) {
                                            finalShared[i][j]++;
                                        }
                                        finalTotal[i][j]++;
                                    }
                                }
                            });
                            processed++;
                        }

                    }
                }
                ln = tf.readLine();
                ctr++;
                if (ctr % 1000 == 0) {
                    System.out.print(ctr + " lines parsed, " + processed + " current vcf: " + vcffile + "\r");
                }
            }
            tf.close();
            System.out.print(ctr + " lines parsed, " + processed + " current vcf: " + vcffile + "\r");

        }

        TextFile out = new TextFile(output, TextFile.W);
        TextFile out2 = new TextFile(output + "-percShared.txt", TextFile.W);
        String header = "-\t" + Strings.concat(samples, Strings.tab);
        out.writeln(header);
        out2.writeln(header);
        for (int i = 0; i < grm.length; i++) {
            for (int j = 0; j < grm[i].length; j++) {
                if (total[i][j] > 0) {
                    grm[i][j] /= total[i][j];
                    grm[j][i] = grm[i][j];
                    shared[i][j] /= total[i][j];
                    shared[j][i] = shared[i][j];
                }
            }
            String ln = samples.get(i) + "\t" + Strings.concat(grm[i], Strings.tab);
            String ln2 = samples.get(i) + "\t" + Strings.concat(shared[i], Strings.tab);
            out.writeln(ln);
            out2.writeln(ln2);
        }
        out.close();
        out2.close();

    }
}
