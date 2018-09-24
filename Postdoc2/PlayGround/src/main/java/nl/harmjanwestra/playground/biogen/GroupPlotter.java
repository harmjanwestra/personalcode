package nl.harmjanwestra.playground.biogen;

import JSci.maths.ArrayMath;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class GroupPlotter {

    public static void main(String[] args) {

        GroupPlotter p = new GroupPlotter();

        String datafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\rnapca\\pc1_2.txt";
        String groupfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-CategoricalMetaData-RNAIDs.txt";
        String outfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\rnapca\\pc1_2.pdf";

        try {
            p.plot(datafile, null, null, groupfile, outfile);
        } catch (IOException e) {
            e.printStackTrace();
        } catch (DocumentException e) {
            e.printStackTrace();
        }

    }

    public void plot(String datafile, Integer col1, Integer col2, String groupfile, String outfile) throws IOException, DocumentException {
        if (col1 == null) {
            col1 = 0;
        }
        if (col2 == null) {
            col2 = 1;
        }
        DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(datafile);


        DoubleMatrix1D col1vals = ds.getCol(col1);
        double minx = ArrayMath.min(col1vals.toArray());
        double maxx = ArrayMath.min(col1vals.toArray());
        DoubleMatrix1D col2vals = ds.getCol(col2);
        double miny = ArrayMath.min(col2vals.toArray());
        double maxy = ArrayMath.min(col2vals.toArray());
        Range range = new Range(minx, miny, maxx, maxy);

        TextFile tf = new TextFile(groupfile, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);

        ArrayList<String> sets = new ArrayList<String>();
        ArrayList<HashMap<String, String>> sampleToSetToGroup = new ArrayList<HashMap<String, String>>();
        ArrayList<HashSet<String>> groupsPerSet = new ArrayList<HashSet<String>>();

        for (int i = 1; i < header.length; i++) {
            sets.add(header[i]);
            sampleToSetToGroup.add(new HashMap<>());
            groupsPerSet.add(new HashSet<>());
        }

        System.out.println(sets.size() + " sets present");
        sets.add("Unlabeled");

        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String sample = elems[0];

            for (int g = 1; g < elems.length; g++) {
                int set = g - 1;
                String group = elems[g];
                sampleToSetToGroup.get(set).put(sample, group);
                groupsPerSet.get(set).add(group);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        sampleToSetToGroup.add(new HashMap<>());
        groupsPerSet.add(new HashSet<>());

        int plotsperrow = 1;
        int nrrows = (int) Math.ceil((double) sets.size() / plotsperrow);

        String unknowngrouplabel = "UnknownGroup";

        for (int set = 0; set < sets.size(); set++) {
            Grid grid = new Grid(500, 500, nrrows, plotsperrow, 100, 100);
            HashSet<String> groupset = groupsPerSet.get(set);
            ArrayList<String> groups = new ArrayList<>();
            groups.addAll(groupset);

            Collections.sort(groups);
            groups.add(unknowngrouplabel);
            HashMap<String, String> sampleToGroup = sampleToSetToGroup.get(set);

            ScatterplotPanel p = new ScatterplotPanel(1, 1);

            // group xy datapoints per sample
            ArrayList<String> samples = ds.getRowObjects();
            // inventorize groups
            HashMap<String, Integer> groupToInt = new HashMap<>();
            ArrayList<ArrayList<Double>> x = new ArrayList<>();
            ArrayList<ArrayList<Double>> y = new ArrayList<>();

            int ctr = 0;
            for (String group : groups) {
                groupToInt.put(group, ctr);
                x.add(new ArrayList<>());
                y.add(new ArrayList<>());
                ctr++;
            }

            for (int s = 0; s < samples.size(); s++) {
                String sample = samples.get(s);
                String group = sampleToGroup.get(sample);

                if (group == null) {
                    group = unknowngrouplabel;
                }

                Integer gindex = groupToInt.get(group);
                x.get(gindex).add(ds.getElementQuick(s, col1));
                y.get(gindex).add(ds.getElementQuick(s, col2));

            }


            // remove groups without data
            ArrayList<ArrayList<Double>> tmpx = new ArrayList<>();
            ArrayList<ArrayList<Double>> tmpy = new ArrayList<>();
            ArrayList<String> tmpgroups = new ArrayList<>();
            for (int g = 0; g < groups.size(); g++) {
                if (!x.get(g).isEmpty()) {
                    tmpx.add(x.get(g));
                    tmpy.add(y.get(g));
                    tmpgroups.add(groups.get(g));
                }
            }


            // convert to double[][];
            double[][] dataX = new double[tmpgroups.size()][];
            double[][] dataY = new double[tmpgroups.size()][];

            for (int g = 0; g < tmpgroups.size(); g++) {
                dataX[g] = Primitives.toPrimitiveArr(tmpx.get(g));
                dataY[g] = Primitives.toPrimitiveArr(tmpy.get(g));
            }

            System.out.println();
            System.out.println(sets.get(set));
            for (int g = 0; g < tmpgroups.size(); g++) {
                System.out.println(tmpgroups.get(g) + "\t" + dataX[g].length + "\t" + dataY[g].length);
            }

            ColorBlindTheme theme = new ColorBlindTheme();
            p.setTheme(theme);
            p.setData(dataX, dataY);
            p.setDatasetLabels(tmpgroups.toArray(new String[0]));
            p.setLabels(ds.getColObjects().get(col1), ds.getColObjects().get(col2));
            p.setTitle(sets.get(set));
//            p.setDataRange(range);
            p.setPlotElems(true, true);
            p.setAlpha(0.4f);
            grid.addPanel(p);
            grid.draw(outfile + "-" + sets.get(set) + ".png");
        }


    }
}
