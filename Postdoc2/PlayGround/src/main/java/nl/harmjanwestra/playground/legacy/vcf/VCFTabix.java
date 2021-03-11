package nl.harmjanwestra.playground.legacy.vcf;


import org.broad.tribble.readers.TabixReader;
import umcg.genetica.features.Feature;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 10/6/16.
 */
public class VCFTabix {

    private final String tabixfile;
    private final TabixReader treader;

    public VCFTabix() {
        this.tabixfile = null;
        this.treader = null;
    }

    public VCFTabix(String filename) throws IOException {
        this.tabixfile = filename;
        if (Gpio.exists(tabixfile + ".tbi")) {
            treader = new TabixReader(tabixfile);
        } else {
            System.out.println("Could not find tabix index: " + tabixfile + ".tbi");
            treader = null;
        }

    }

    public boolean[] getSampleFilter(String tabixsamplelimit) throws IOException {
        return getSampleFilter(tabixsamplelimit, this.tabixfile);
    }

    public boolean[] getSampleFilter(String tabixsamplelimit, String vcffile) throws IOException {
        if (tabixsamplelimit == null) {
            return null;
        }
        boolean[] samplesToInclude = null;
        VCFGenotypeData d = new VCFGenotypeData(vcffile);
        ArrayList<String> tabixSamples = d.getSamples();
        d.close();

        samplesToInclude = new boolean[tabixSamples.size()];

        TextFile tf = new TextFile(tabixsamplelimit, TextFile.R);
        String[] elems = tf.readLineElems(Strings.whitespace);
        HashSet<String> allSamplesInListSet = new HashSet<String>();
        while (elems != null) {
            allSamplesInListSet.add(elems[0]);
            elems = tf.readLineElems(Strings.whitespace);
        }
        tf.close();
        for (int s = 0; s < tabixSamples.size(); s++) {
            if (allSamplesInListSet.contains(tabixSamples.get(s))) {
                samplesToInclude[s] = true;
            }
        }
        return samplesToInclude;
    }

    public TabixReader.Iterator query(Feature region) throws IOException {
        int start = region.getStart();
        if (start < 1) {
            start = 1;
        }
        int stop = region.getStop();

//		System.out.println("Query: " + region.getChromosome().getNumber() + ":" + start + "-" + stop + "\t" + treader.getSource());
        String queryStr = region.getChromosome().getNumber() + ":" + start + "-" + stop;
        try {

            TabixReader.Iterator window = treader.query(queryStr);
            return window;
        } catch (ArrayIndexOutOfBoundsException e) {
            queryStr = "chr" + region.getChromosome().getNumber() + ":" + start + "-" + stop;
            TabixReader.Iterator window = treader.query(queryStr);
            return window;
        }
    }

    public void close() {
        treader.close();
    }

    public ArrayList<VCFVariant> getAllVariants(Feature f, boolean[] samplefilter) throws IOException {
        TabixReader.Iterator iterator = query(f);
        String next = iterator.next();
        ArrayList<VCFVariant> output = new ArrayList<>();
        while (next != null) {
            output.add(new VCFVariant(next, VCFVariant.PARSE.ALL, samplefilter));
            next = iterator.next();
        }
        return output;
    }

    public VCFVariant getVariant(Feature snpFeature, boolean[] samplefilter) throws IOException {
        Feature s = new Feature(snpFeature);
        s.setStart(s.getStart() - 1);
        s.setStop(s.getStop() + 1);


        TabixReader.Iterator it = query(s);
//        if (snpFeature.getStart() == 82082973 || snpFeature.getStart() == 81359709) {
//            System.out.println(snpFeature + " --> " + s + "\twindow found: " + (it != null));
//        }
        if (it != null) {

            String str = it.next();

            while (str != null) {
                VCFVariant variant = new VCFVariant(str, VCFVariant.PARSE.HEADER, samplefilter);
//                if (snpFeature.getStart() == 82082973 || snpFeature.getStart() == 81359709) {
//                    System.out.println(variant.asFeature());
//                }
                if (variant.asFeature().overlaps(snpFeature)) {
                    variant = new VCFVariant(str, VCFVariant.PARSE.ALL, samplefilter);
//                    if (snpFeature.getStart() == 82082973 || snpFeature.getStart() == 81359709) {
//                        System.out.println(variant.asFeature() + " found!");
//                    }
                    return variant;
                }
                str = it.next();
            }
        }

        return null;
    }
}
