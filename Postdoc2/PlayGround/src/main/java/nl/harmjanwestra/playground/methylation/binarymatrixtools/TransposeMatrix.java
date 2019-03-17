package nl.harmjanwestra.playground.methylation.binarymatrixtools;

import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;

public class TransposeMatrix {

    public static void main(String[] args) {
        String file = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\test.txt";
        String out = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\test.dat";
        String out2 = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\test2.txt";

        MatrixConverter cv = new MatrixConverter();
        TransposeMatrix t = new TransposeMatrix();
        String outtranspose = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\test-transpose.dat";
        try {
            t.transposeLargeMatrix(out, outtranspose, 2);
            String outtransposetext = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\test-transpose.txt";
            cv.toText(outtranspose, outtransposetext);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    public void transpose(String in, String out) throws IOException {
        System.out.println("Transpose.");
        System.out.println(in);
        System.out.println(out);
        // Gpio.delete(out);
        DiskBasedBinaryDoubleMatrixReader reader = new DiskBasedBinaryDoubleMatrixReader(new File(in));
        DiskBasedBinaryDoubleMatrixWriter writer = new DiskBasedBinaryDoubleMatrixWriter();
        writer.initializeFullMatrix(reader.colIds, reader.rowIds, out);
        writer.close();
        writer.open(out);
        for (int i = 0; i < reader.nrRows; i++) {
            double[] row = reader.getNextRow();
            System.out.println(Strings.concat(row, Strings.tab));
            for (int j = 0; j < reader.nrCols; j++) {

                writer.write(j, i, row[j]);
            }
        }
        writer.close();
    }

    public void transposeLargeMatrix(String in, String out, int rowBufferSize) throws IOException {
        System.out.println("Transposing big matrix.");
        System.out.println(in);
        System.out.println(out);

        // Gpio.delete(out);
        DiskBasedBinaryDoubleMatrixReader reader = new DiskBasedBinaryDoubleMatrixReader(new File(in));
        DiskBasedBinaryDoubleMatrixWriter writer = new DiskBasedBinaryDoubleMatrixWriter();
        writer.initializeFullMatrix(reader.colIds, reader.rowIds, out);
        writer.close();
        writer.open(out);


        System.out.println("Stuff is printed here");
        double[][] buffer = new double[rowBufferSize][];
        double[][] transpose = new double[reader.colIds.length][rowBufferSize];

        int ctr = 0;

        int buffernr = 0;
        for (int i = 0; i < reader.nrRows; i++) {
            if (ctr == rowBufferSize) {
                // buffer full
                for (int q = 0; q < buffer.length; q++) {
                    for (int z = 0; z < buffer[q].length; z++) {
                        transpose[z][q] = buffer[q][z];
                    }
                }

                // data has all columns
                int colstart = buffernr * rowBufferSize;
                for (int q = 0; q < transpose.length; q++) {
                    writer.writeBlock(q, colstart, transpose[q]);
                }

                buffernr++;
                ctr = 0;
            }
            if (ctr != rowBufferSize) {
                // space left in the buffer.
                double[] row = reader.getNextRow();
                buffer[ctr] = row;
                ctr++;

            }


        }

        if (ctr > 0) {
            transpose = new double[reader.colIds.length][ctr];

            // buffer full
            for (int q = 0; q < ctr; q++) {
                for (int z = 0; z < buffer[q].length; z++) {

                    transpose[z][q] = buffer[q][z];
                }
            }

            // data has all columns,
            int colstart = buffernr * rowBufferSize;
            for (int q = 0; q < transpose.length; q++) {
                writer.writeBlock(q, colstart, transpose[q]);
            }
        }
        writer.close();


    }

}
