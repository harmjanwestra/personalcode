package nl.harmjanwestra.playground.methylation;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRandomAccessReader;

import java.io.IOException;

public class CompareBinaryMatrices {


    public void run(String input1, String input2) throws IOException {


        DoubleMatrixDatasetRandomAccessReader reader1 = new DoubleMatrixDatasetRandomAccessReader(input1);
        DoubleMatrixDatasetRandomAccessReader reader2 = new DoubleMatrixDatasetRandomAccessReader(input2);

        boolean ok = true;
        if (reader1.rows() != reader2.rows()) {
            ok = false;
            System.out.println(input1 + " and " + input2 + " have different numbers of rows");

        }
        if (reader1.cols() != reader2.cols()) {
            ok = false;
            System.out.println(input1 + " and " + input2 + " have different numbers of cols");
        }
        if (ok) {
            ProgressBar pb = new ProgressBar(reader1.rows());
            for (int r = 0; r < reader1.rows(); r++) {
                double[] row1 = reader1.getRow(r);
                double[] row2 = reader2.getRow(r);

                // compare rows
                double diffsum = 0;
                for (int c = 0; c < row1.length; c++) {
                    diffsum += (row1[c] - row2[c]);
                }
                if (diffsum > 0.0001 || diffsum < -0.0001) {
                    System.out.println("Row " + r + " difference: " + diffsum);
                }
                pb.iterate();
            }
            pb.close();
        }

    }
}
