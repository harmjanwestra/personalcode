package nl.harmjanwestra.playground.methylation;

import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
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

    }

    public void run(String projectfile, String indir, String outloc) throws IOException {

        HashSet<String> projects = new HashSet<>();
        TextFile tf = new TextFile(projectfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String project = elems[0];
            projects.add(project);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        for (String project : projects) {
            System.out.println("Parsing project: " + project);
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
                        String probe = Strings.cache(elems[0]);
                        if (!probemap.containsKey(probe)) {
                            probemap.put(probe, probes.size());
                            probes.add(probe);
                        }

                        Double val = Double.parseDouble(elems[1]);

                        Pair<String, Double> d = new Pair<String, Double>(probe, val);
                        indData.add(d);
                        melems = tf2.readLineElems(TextFile.tab);
                    }
                    tf2.close();
                    data.add(indData);
                    inds.add(f.getName());
                }
            }

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

                    out.save(outloc + project + "-beta.txt.gz");

                    // convert to M-values
                    double[][] raw = out.getMatrix().toArray();
                    ConvertBetaAndMvalues.transformToMvalue(raw);
                    out.setMatrix(raw);
                    out.save(outloc + project + "-Mval.txt.gz");
                } catch (Exception e) {
                    e.printStackTrace();
                }


            }


        }


    }

}
