package nl.harmjanwestra.playground.methylation;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.methylation.ConvertBetaAndMvalues;
import umcg.genetica.methylation.ConverteUandMTo;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class ProcessTCGA {

    public static void main(String[] args) {

        String dir = "S:\\projects\\2018-tcga\\solidtissue\\dnameth\\";
        String projectfile = "S:\\projects\\2018-tcga\\inputmeth\\gdc_sample_sheet.2018-08-06.tsv";
        String outloc = "S:\\projects\\2018-tcga\\solidtissue\\dnamethmerged\\";

        ProcessTCGA t = new ProcessTCGA();
        try {
//            t.runMeth(projectfile, dir, outloc);

            dir = "S:\\projects\\2018-tcga\\solidtissue\\rna\\";
            projectfile = "S:\\projects\\2018-tcga\\inputrna\\gdc_sample_sheet.2018-08-06.tsv";
            outloc = "S:\\projects\\2018-tcga\\solidtissue\\rnamerged\\";
            t.run(projectfile, dir, outloc);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void run(String projectfile, String indir, String outloc) throws IOException {

        HashSet<String> projects = new HashSet<>();
        TextFile tf = new TextFile(projectfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            if (elems[7].contains("Solid")) {
                String project = elems[4];
                projects.add(project);
            }

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        int pctr = 0;
        for (String project : projects) {
            System.out.println("Parsing project: " + project + "\t" + pctr + "/" + projects.size());
            String folder = indir + project;
            File[] files = new File(folder).listFiles();
            ArrayList<File> selected = new ArrayList<>();
            HashMap<String, Integer> probemap = new HashMap<>();
            ArrayList<String> probes = new ArrayList<>();

            ArrayList<ArrayList<Pair<String, Double>>> data = new ArrayList<>();
            ArrayList<String> inds = new ArrayList<>();
            for (File f : files) {
                if (f.isFile() && f.getName().endsWith(".gz")) {

                    selected.add(f);

                    ArrayList<Pair<String, Double>> indData = new ArrayList<>();
                    TextFile tf2 = new TextFile(f, TextFile.R);
                    tf2.readLine();
                    String[] melems = tf2.readLineElems(TextFile.tab);
                    while (melems != null) {
                        String probe = Strings.cache(melems[0]);
                        if (!probemap.containsKey(probe)) {
                            probemap.put(probe, probes.size());
                            probes.add(probe);
                        }

                        Double val = Double.NaN;
                        if (!melems[1].contains("NA")) {
                            val = Double.parseDouble(melems[1]);
                        }

                        Pair<String, Double> d = new Pair<String, Double>(probe, val);
                        indData.add(d);
                        melems = tf2.readLineElems(TextFile.tab);
                    }
                    String indname = f.getName().replaceAll(".txt.gz", "");
                    System.out.println(project + "\t" + inds.size() + "/" + files.length + "\t" + indname + "\t" + indData.size() + " probes");
                    tf2.close();
                    data.add(indData);
                    inds.add(indname);
                }

            }

            System.out.println(project + "\t" + inds.size() + " individuals");

            if (!inds.isEmpty()) {
                DoubleMatrixDataset<String, String> out = new DoubleMatrixDataset<>(probes.size(), inds.size());
                try {
                    out.setRowObjects(probes);
                    out.setColObjects(inds);

                    for (int i = 0; i < data.size(); i++) {
                        ArrayList<Pair<String, Double>> d = data.get(i);
                        for (Pair<String, Double> p : d) {
                            Integer pid = probemap.get(p.getLeft());
                            out.getMatrix().setQuick(pid, i, p.getRight());
                        }
                    }
                    System.out.println("Saving: " + outloc + project + "-beta.txt.gz");
                    out.save(outloc + project + "-beta.txt.gz");

                    // convert to M-values
                    double[][] raw = out.getMatrix().toArray();
                    ConvertBetaAndMvalues.transformToMvalue(raw);
                    out.setMatrix(raw);
                    System.out.println("Saving: " + outloc + project + "-Mval.txt.gz");
                    out.save(outloc + project + "-Mval.txt.gz");

                } catch (Exception e) {
                    e.printStackTrace();
                }


            }
            pctr++;

        }


    }

}
