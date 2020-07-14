package nl.harmjanwestra.methylation.binarymatrixtools;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixConverter;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;

public class TransposeMatrix {

    public static void main(String[] args) {

        try {
            String matrix = "D:\\tmp\\matrixtest\\testmatrix.txt";
            String matrixbinary = "D:\\tmp\\matrixtest\\testmatrix-binary";
            TransposeMatrix t = new TransposeMatrix();
            t.createMatrix(matrix, 100);
            DoubleMatrixConverter.TextToBinary(matrix, matrixbinary);

            String matrixbinarytp = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed";
            String matrixbinarytptxt = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed.txt";

            t.transposeLargeMatrix(matrixbinary, matrixbinarytp, 10);
            DoubleMatrixConverter.BinaryToText(matrixbinarytp + ".dat", matrixbinarytptxt);

            String matrixbinarytptp = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed-transposed";
            t.transposeLargeMatrix(matrixbinarytp, matrixbinarytptp, 10);
            String matrixbinarytptptxt = "D:\\tmp\\matrixtest\\testmatrix-binary-transposed-transposed.txt";

            DoubleMatrixConverter.BinaryToText(matrixbinarytptp + ".dat", matrixbinarytptptxt);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void createMatrix(String output, int size) throws IOException {
        double[][] matrix = new double[size][size];
        int ctr = 0;
        String[] rows = new String[size];
        String[] cols = new String[size];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                matrix[i][j] = ctr;
                ctr++;
            }
            rows[i] = "Row-" + i;
            cols[i] = "Col-" + i;
        }

        TextFile tf = new TextFile(output, TextFile.W);
        String header = "-\t" + Strings.concat(cols, Strings.tab);
        tf.writeln(header);
        for (int r = 0; r < size; r++) {
            tf.writeln("Row" + r + "\t" + Strings.concat(matrix[r], Strings.tab));
        }
        tf.close();

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


//        System.out.println("Stuff is printed here");
        double[][] buffer = new double[rowBufferSize][];
        double[][] transpose = new double[reader.colIds.length][rowBufferSize];

        int ctr = 0;

        int buffernr = 0;
        ProgressBar pb = new ProgressBar(reader.nrRows, "Transposing...");
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

            pb.iterate();

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

        pb.close();

    }

}
