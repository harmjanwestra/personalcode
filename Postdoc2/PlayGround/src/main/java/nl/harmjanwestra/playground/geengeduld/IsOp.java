package nl.harmjanwestra.playground.geengeduld;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.text.Strings;

import java.util.ArrayList;
import java.util.stream.IntStream;

public class IsOp {

    public static void main(String[] args) {
        IsOp c = new IsOp();
        if (args.length < 2) {
            System.out.println("Usage: input.txt.gz output.txt.gz");
            System.exit(-1);
        }
        try {
            c.run(args[0], args[1]);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public void run(String dsfile, String outfile) throws Exception {
        System.out.println("Loading: " + dsfile);
        DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleTextData(dsfile, '\t');


        double[][] matrix = new double[ds.columns()][ds.columns()];
        ProgressBar pb = new ProgressBar(ds.columns(), "Correlating! ");

        IntStream.range(0, ds.columns()).parallel().forEach(c -> {
            double[] x = ds.viewCol(c).toArray();
            for (int c2 = c + 1; c2 < ds.columns(); c2++) {
                double[] y = ds.viewCol(c2).toArray();
                Pair<Double, Integer> result = Correlation.correlateWithNaNValues(x, y);
                matrix[c][c2] = result.getLeft();
            }
            pb.iterateSynched();
        });
        pb.close();

        for (int i = 0; i < matrix.length; i++) {
            matrix[i][i] = 1;
            for (int j = i + 1; j < matrix.length; j++) {
                matrix[j][i] = matrix[i][j];
            }
        }

        pb.close();

        System.out.println("Writing: " + outfile);

        ArrayList<String> ids = ds.getColObjects();

        TextFile tf = new TextFile(outfile, TextFile.W);
        String header = "-\t" + Strings.concat(ids, Strings.tab);
        tf.writeln(header);
        for (int r = 0; r < matrix.length; r++) {
            String outln = ids.get(r) + "\t" + Strings.concat(matrix[r], Strings.tab);
            tf.writeln(outln);
        }
        tf.close();

    }


}
