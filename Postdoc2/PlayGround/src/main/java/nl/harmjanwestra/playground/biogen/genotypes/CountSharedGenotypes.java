package nl.harmjanwestra.playground.biogen.genotypes;

import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class CountSharedGenotypes {

    public static void main(String[] args) {

    }

    public void run(String input, String output) throws IOException {

        HashMap<String, Integer> sampleHash = new HashMap<>();
        ArrayList<String> samples = new ArrayList<>();

        double[][] shared = null;
        int[][] total = null;

        for (int chr = 1; chr < 23; chr++) {
            TextFile tf = new TextFile(input, TextFile.R);
            String ln = tf.readLine();
            int[] sampleOrder = null;
            while (ln != null) {
                if (ln.startsWith("#")) {
                    // skip
                } else if (ln.startsWith("#CHROM")) {
                    // header line

                    // initialize if not yet done so
                    if (total == null) {
                        String[] elems = ln.split("\t");
                        for (int i = 9; i < elems.length; i++) {
                            sampleHash.put(elems[i], i - 9);
                            samples.add(elems[i]);
                        }
                        sampleOrder = new int[samples.size()];
                        for (int i = 0; i < samples.size(); i++) {
                            sampleOrder[i] = i;
                            total = new int[i][samples.size() - i]; // create staggered array of arrays
                            shared = new double[i][samples.size() - i];
                        }
                    } else {
                        // reorder samples to first file
                        String[] elems = ln.split("\t");
                        sampleOrder = new int[samples.size()];
                        Arrays.fill(sampleOrder, -1);
                        for (int i = 9; i < elems.length; i++) {
                            Integer id = sampleHash.get(elems[i]);
                            if (id != null) {
                                sampleOrder[i] = id;
                            }
                        }
                    }
                } else {
                    // variant line.
                    // parse variant
                    VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
                    byte[] dosages = variant.getGenotypesAsByteVector();

                }
                ln = tf.readLine();
            }
            tf.close();
        }


    }
}
