package nl.harmjanwestra.methylation.binarymatrixtools;


import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class MatrixConverter {

    public static void main(String[] args) {
        String filein = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\GPL13534_MergedFromIdat.txt.gz";
        String fileout = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\GPL13534_MergedFromIdat.dat";
        MatrixConverter c = new MatrixConverter();
        try {
            c.toBinary(filein, fileout);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void toBinary(String in, String out) throws IOException {
        System.out.println("Converting matrix to binary..");
        System.out.println("In " + in);
        System.out.println("Out: " + out);
        BinaryFile bf = new BinaryFile(out, BinaryFile.W);

        ArrayList<String> rowIds = new ArrayList<String>();
        System.out.println("Gathering rows in file: " + in);
        int ctr2 = 0;
        TextFile tf = new TextFile(in, TextFile.R);
        tf.readLine();
        String ln = tf.readLine();
        while (ln != null) {
            if (ln.trim().length() > 0) {
                String[] elems = Strings.subsplit(ln, Strings.tab, 0, 1);
                rowIds.add(elems[0]);
            }
            ln = tf.readLine();
            ctr2++;
            if (ctr2 % 10 == 0) {
                System.out.print("\r" + ctr2 + " lines parsed");
            }
        }
        tf.close();
        tf.open();
        System.out.println();
        System.out.println(rowIds.size() + " total rows");
        bf.writeInt(rowIds.size());
        for (int i = 0; i < rowIds.size(); i++) {
            bf.writeString(rowIds.get(i));
        }

        String[] header = tf.readLineElems(TextFile.tab);
        System.out.println((header.length - 1) + " columns.");
        bf.writeInt(header.length - 1);
        for (int i = 1; i < header.length; i++) {
            bf.writeString(header[i]);
        }


        String[] elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;
        ProgressBar pb = new ProgressBar(rowIds.size(), "Writing binary doubles");
        while (elems != null) {
            if (elems[0].trim().length() > 0) {
                for (int i = 1; i < elems.length; i++) {
                    bf.writeDouble(Double.parseDouble(elems[i]));
                }
                pb.iterate();

            }
            elems = tf.readLineElems(TextFile.tab);
        }
        pb.close();
        tf.close();
        bf.close();
    }

    public void toText(String input, String output) throws IOException {
        System.out.println("Converting binary to text..");
        System.out.println("In " + input);
        System.out.println("Out: " + output);
        TextFile tf = new TextFile(output, TextFile.W);
        DiskBasedBinaryDoubleMatrixReader bf = new DiskBasedBinaryDoubleMatrixReader(new File(input));

        String header = "-\t" + Strings.concat(bf.getCols(), Strings.tab);
        tf.writeln(header);
        double[] buffer = new double[bf.getNrCols()];
        String[] rowIds = bf.getRowIds();
        int nrRows = bf.getNrRows();
        for (int row = 0; row < nrRows; row++) {
            String ln = rowIds[row] + "\t" + Strings.concat(
                    bf.getRow(row, buffer), Strings.tab);
            if (ln.trim().length() > 0) {
                tf.writeln(ln);
            }
        }
        bf.close();
        tf.close();

    }
}
