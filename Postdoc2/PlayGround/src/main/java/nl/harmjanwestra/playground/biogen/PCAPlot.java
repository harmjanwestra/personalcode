package nl.harmjanwestra.playground.biogen;

import com.itextpdf.text.DocumentException;
import umcg.genetica.containers.Triple;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.graphics.themes.DefaultTheme;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class PCAPlot {

    public static void main(String[] args) {
//        String knownClasses = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p3v5a-sampleinfo.txt";
//        String superclasses = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p3v5a-superpopulations.txt";


//		String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\genotypepca\\plink.eigenvec.txt";
//		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\genotypepca\\AMPAD-";
//		String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\genotypepca\\plink.eigenvec";
//		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\genotypepca\\GTEx-";
//        String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\plink.eigenvec.txt";
//        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-";
//        String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\genotypepca\\plink.eigenvec.txt";
//        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\genotypepca\\CMC-";
//        String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\dnapca\\plink.eigenvec.txt";
//        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\dnapca\\Braineac-";


        PCAPlot p = new PCAPlot();

        try {
//            String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\hap550qc\\lng_4tissue_brain_qtl_fix-pca - Copy.eigenvec.txt";
//            String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\hap550-";
//
//            p.plot(pcafile, knownClasses, superclasses, output);
//            pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\hap660qc\\eQTL_v2_fix-pca.eigenvec";
//            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\hap610-";
//
//            p.plot(pcafile, knownClasses, superclasses, output);

            String knownClasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-tissue.annot.txt";
            String superclasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-tissue.annot-groups.txt";
            String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\mergednorm\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.PCAOverSamplesEigenvectors.txt.gz";
//            String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\mergednorm\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors.txt.gz";
            String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\MergedRNASeqPCAPlot-Centered-PC1v2-tissue.pdf";
            p.plot(pcafile, knownClasses, superclasses, 1, 2, 7, output);
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\MergedRNASeqPCAPlot-Centered-PC3v4-tissue.pdf";
            p.plot(pcafile, knownClasses, superclasses, 3, 4, 7, output);

            knownClasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-dataset.annot.txt";
            superclasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-dataset.annot-groups.txt";
            pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\mergednorm\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.PCAOverSamplesEigenvectors.txt.gz";
//            String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\mergednorm\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors.txt.gz";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\MergedRNASeqPCAPlot-Centered-PC1v2-dataset.pdf";
            p.plot(pcafile, knownClasses, superclasses, 1, 2, 7, output);
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\MergedRNASeqPCAPlot-Centered-PC3v4-dataset.pdf";
            p.plot(pcafile, knownClasses, superclasses, 3, 4, 7, output);


        } catch (IOException e) {
            e.printStackTrace();
        } catch (DocumentException e) {
            e.printStackTrace();
        }

    }


    public void plot(String pcafile, String knownclasses, String superclasses, int col1, int col2, int startk, String output) throws IOException, DocumentException {

        HashMap<String, String> superpop = new HashMap<>();
        TextFile tfc = new TextFile(superclasses, TextFile.R);
        tfc.readLine();
        String[] elemsf = tfc.readLineElems(TextFile.tab);
        while (elemsf != null) {
            superpop.put(elemsf[0], elemsf[2]);
            elemsf = tfc.readLineElems(TextFile.tab);
        }
        tfc.close();

        String nullpop = "NotDefined";


        HashMap<String, HashSet<String>> samplesPerPop = new HashMap<>();
        HashMap<String, String> sampleToPop = new HashMap<>();
        TextFile tf = new TextFile(knownclasses, TextFile.R);

        tf.readLine();

        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            String sample = elems[0];
            String type = elems[1];
            type = superpop.get(type);
            HashSet<String> set = samplesPerPop.get(type);
            if (set == null) {
                set = new HashSet<>();
            }
            set.add(sample);
            samplesPerPop.put(type, set);
            sampleToPop.put(sample, type);

            elems = tf.readLineElems(Strings.whitespace);

        }
        tf.close();

        System.out.println(sampleToPop.size() + " with known population.");
        System.out.println(samplesPerPop.size() + " populations total.");


        ArrayList<String> populations = new ArrayList<>();
        ArrayList<String> populationstmp = new ArrayList<>();
        populationstmp.addAll(samplesPerPop.keySet());
        Collections.sort(populationstmp);
        populations.add(nullpop);
        populations.addAll(populationstmp);


        HashMap<String, Integer> populationIndex = new HashMap<>();
        ArrayList<ArrayList<Double>> xvals = new ArrayList<>();
        ArrayList<ArrayList<Double>> yvals = new ArrayList<>();
        ArrayList<ArrayList<Double>> xvalsunlabeled = new ArrayList<>();
        ArrayList<ArrayList<Double>> yvalsunlabeled = new ArrayList<>();
        for (int i = 0; i < populations.size(); i++) {
            String pop = populations.get(i);
            populationIndex.put(pop, i);
            xvals.add(new ArrayList<>());
            xvalsunlabeled.add(new ArrayList<>());
            yvals.add(new ArrayList<>());
            yvalsunlabeled.add(new ArrayList<>());
        }

        ArrayList<Triple<String, Double, Double>> data = new ArrayList<Triple<String, Double, Double>>();
        TextFile tf2 = new TextFile(pcafile, TextFile.R);

        String ln = tf2.readLine();
        double minX = Double.MAX_VALUE;
        double maxX = -Double.MAX_VALUE;
        double minY = Double.MAX_VALUE;
        double maxY = -Double.MAX_VALUE;
        int ctr = 0;
        while (ln != null) {
            String pcadata = ln;
            while (pcadata.contains("  ")) {
                ln.replaceAll("  ", " ");
            }
            String[] pcaelems = Strings.whitespace.split(pcadata);
            String sample = pcaelems[0];
            try {
                Double pca1 = Double.parseDouble(pcaelems[col1]);

                Double pca2 = Double.parseDouble(pcaelems[col2]);
                if (!Double.isNaN(pca1)) {
                    if (pca1 > maxX) {
                        maxX = pca1;
                    }
                    if (pca1 < minX) {
                        minX = pca1;
                    }
                }
                if (!Double.isNaN(pca2)) {
                    if (pca2 > maxY) {
                        maxY = pca2;
                    }
                    if (pca2 < minY) {
                        minY = pca2;
                    }
                }

                data.add(new Triple<>(sample, pca1, pca2));

                String population = sampleToPop.get(sample);
                if (population == null) {
                    population = nullpop;
                }
                sampleToPop.put(sample, population);
                Integer popindex = populationIndex.get(population);
                if (population.equals(nullpop)) {
                    xvalsunlabeled.get(popindex).add(pca1);
                    yvalsunlabeled.get(popindex).add(pca2);
                } else {
                    xvals.get(popindex).add(pca1);
                    yvals.get(popindex).add(pca2);
                }
            } catch (NumberFormatException e) {
                System.out.println("Error parsing line: " + e.getMessage());
                System.out.println("Probably a header I didnt expect.");
            }
            ln = tf2.readLine();
            ctr++;
            if (ctr % 1000 == 0) {
                System.out.println(ctr + " lines parsed");
            }
        }
        tf2.close();

        DefaultTheme def = new DefaultTheme();
        Range r = new Range(minX, minY, maxX, maxY);

        Color[] colors = new Color[]{
                new Color(0, 0, 0),
                new Color(145, 145, 145),

                Color.decode("#f44336".toUpperCase()),
                Color.decode("#2196f3".toUpperCase()),
                Color.decode("#ffeb3b".toUpperCase()),
                Color.decode("#e91e63".toUpperCase()),
                Color.decode("#009688".toUpperCase()),
                Color.decode("#8bc34a".toUpperCase()),
                Color.decode("#ff5722".toUpperCase()),

                Color.decode("#3f51b5".toUpperCase()),
                Color.decode("#00bcd4".toUpperCase()),
                Color.decode("#9c27b0".toUpperCase()),
                Color.decode("#673ab7".toUpperCase()),

                Color.decode("#03a9f4".toUpperCase()),
                Color.decode("#607d8b".toUpperCase()),
                Color.decode("#4caf50".toUpperCase()),
                Color.decode("#cddc39".toUpperCase()),
                Color.decode("#ffc107".toUpperCase()),
                Color.decode("#ff9800".toUpperCase()),
                Color.decode("#795548".toUpperCase()),

        };
        def.setColors(colors);

        // scatterplot
        ScatterplotPanel p1 = preparePanel(r, populations, xvals, yvals, "Reference", "PC" + col1, "PC" + col2, def);
        ScatterplotPanel p2 = preparePanel(r, populations, xvalsunlabeled, yvalsunlabeled, "Dataset", "PC" + col1, "PC" + col2, def);

        Grid grid = new Grid(500, 500, 2, 3, 100, 100);
        grid.addPanel(p1);
        grid.addPanel(p2);

        // knn population assignment
        HashMap<String, String> sampleToPopTmp = assignPopulation(data, sampleToPop, nullpop, startk);

        xvalsunlabeled = new ArrayList<>();
        yvalsunlabeled = new ArrayList<>();
        for (int i = 0; i < populations.size(); i++) {
            xvalsunlabeled.add(new ArrayList<>());
            yvalsunlabeled.add(new ArrayList<>());
        }

        TextFile assignmentout = new TextFile(output + "SampleAssignment.txt", TextFile.W);

        for (Triple<String, Double, Double> sample : data) {
            String prevpop = sampleToPop.get(sample.getLeft());
            if (prevpop.equals(nullpop)) {
                String pop = sampleToPopTmp.get(sample.getLeft());
                Integer popindex = populationIndex.get(pop);
                if (popindex == null) {
                    System.out.println("Could not find population : " + pop + " for sample: " + sample.getLeft());
                } else {
                    xvalsunlabeled.get(popindex).add(sample.getMiddle());
                    yvalsunlabeled.get(popindex).add(sample.getRight());
                    assignmentout.writeln(sample.getLeft() + "\t" + pop);
                }
            }
        }
//        for (Triple<String, Double, Double> sample : data2) {
//            String prevpop = sampleToPop.get(sample.getLeft());
//            if (prevpop.equals(nullpop)) {
//                String pop = sampleToPopTmp.get(sample.getLeft());
//                Integer popindex = populationIndex.get(pop);
//                if (popindex == null) {
//                    System.out.println("Could not find population : " + pop + " for sample: " + sample.getLeft());
//                } else {
//                    assignmentout.writeln(sample.getLeft() + "\t" + pop);
//                }
//            }
//        }

        assignmentout.close();

        ScatterplotPanel p3 = preparePanel(r, populations, xvalsunlabeled, yvalsunlabeled, "Dataset", "PC" + col1, "PC" + col2, def);
        grid.addPanel(p3);


        grid.draw(output + "SampleAssignment.pdf");
    }

    private ScatterplotPanel preparePanel(Range datarange,
                                          ArrayList<String> populations,
                                          ArrayList<ArrayList<Double>> x,
                                          ArrayList<ArrayList<Double>> y,
                                          String title,
                                          String xlabel,
                                          String ylabel,
                                          DefaultTheme def) {
        double[][] xprim = new double[populations.size()][];
        double[][] yprim = new double[populations.size()][];
        for (int i = 0; i < x.size(); i++) {
            xprim[i] = Primitives.toPrimitiveArr(x.get(i));
            yprim[i] = Primitives.toPrimitiveArr(y.get(i));
        }

        ScatterplotPanel punlab = new ScatterplotPanel(1, 1);
        punlab.setDataRange(datarange);
        punlab.setData(xprim, yprim);
        punlab.setTitle(title);
        punlab.setDatasetLabels(populations.toArray(new String[0]));
        punlab.setLabels(xlabel, ylabel);
        punlab.setPlotElems(true, true);
        punlab.setTheme(def);
        punlab.setAlpha(0.8f);

        return punlab;
    }

    private HashMap<String, String> assignPopulation(ArrayList<Triple<String, Double, Double>> data, HashMap<String, String> sampleToPop, String nullpop, int startk) {
        int samplesWithoutPop = 0;
        HashMap<String, String> sampleToPopTmp = new HashMap<>();

        for (Triple<String, Double, Double> sample : data) {
            String pop = sampleToPop.get(sample.getLeft());
            if (pop.equals(nullpop)) {
                int k = startk;
                boolean resolved = false;
                while (!resolved) {
                    String[] neighbors = new String[k];
                    double[] distances = new double[k];
                    double maxdist = -Double.MAX_VALUE;
                    int n = 0;


                    // this sample needs assignment
                    // get startk nearest neighbors that do have an assignment
                    for (Triple<String, Double, Double> sample2 : data) {
                        String pop2 = sampleToPop.get(sample2.getLeft());
                        if (!pop2.equals(nullpop)) {
                            double distance = Math.sqrt(sq(sample.getMiddle() - sample2.getMiddle()) + sq(sample.getRight() - sample2.getRight()));

                            if (n < k) {
                                neighbors[n] = sample2.getLeft();
                                distances[n] = distance;
                                if (distance > maxdist) {
                                    maxdist = distance;
                                }
                                n++;
                            } else if (distance < maxdist) {
                                boolean update = false;
                                double newmax = -Double.MAX_VALUE;

                                for (int i = 0; i < neighbors.length; i++) {
                                    if (!update && (distances[i] > distance)) {
                                        distances[i] = distance;
                                        neighbors[i] = sample2.getLeft();
                                        update = true;
                                    }

                                    if (distances[i] > newmax) {
                                        newmax = distances[i];
                                    }
                                }
                                maxdist = newmax;
                            }
                        }
                    }

                    HashMap<String, Integer> ctr = new HashMap<>();
                    HashMap<String, Double> sumd = new HashMap<>();

                    for (int i = 0; i < neighbors.length; i++) {
                        String pop2 = sampleToPop.get(neighbors[i]);
                        Integer ct = ctr.get(pop2);
                        Double sum = sumd.get(pop2);
                        if (ct == null) {
                            ct = 1;
                            sum = distances[i];
                        } else {
                            ct++;
                            sum += distances[i];
                        }
                        ctr.put(pop2, ct);
                        sumd.put(pop2, sum);
                    }


                    int maxct = 0;
                    String maxpop = null;

                    for (String key : ctr.keySet()) {
                        Integer ct = ctr.get(key);
                        if (ct > maxct) {
                            maxpop = key;
                            maxct = ct;
                        } else if (ct == maxct) {
                            // break ties using absolute sum of distance
                            double sum = sumd.get(key);
                            if (sum < sumd.get(maxpop)) {
                                maxpop = key;
                                maxct = ct;
                            }
                        }
                    }

                    // assign
                    sampleToPopTmp.put(sample.getLeft(), maxpop);
                    resolved = true;
                }
            } else {

                sampleToPopTmp.put(sample.getLeft(), pop);
            }
        }

        return sampleToPopTmp;
    }

    private double sq(double v) {
        return v * v;
    }


}
