package nl.harmjanwestra.playground.methylation;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class MatrixSubset {

    public static void main(String[] args) {
        MatrixSubset s = new MatrixSubset();
        String inmat = args[0]; // "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged\\FromIdat.txt.gz";
        String outmat = args[1]; // "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged\\FromIdat-eQTMProbes.txt.gz";
        String probes = args[2]; // "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2018-05-15-RNACpGCorrelationComparison\\eqtmprobes.txt";
        try {
            s.run(inmat, outmat, probes);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String inmat, String outmat, String probesToKeep) throws IOException {
        TextFile tf = new TextFile(probesToKeep, TextFile.R);
        String ln = tf.readLine();
        HashSet<String> keep = new HashSet<String>();
        while (ln != null) {
            keep.add(ln);
            ln = tf.readLine();
        }
        tf.close();
        System.out.println(keep.size() + " unique probes to keep...");

        TextFile matrix = new TextFile(inmat, TextFile.R);
        TextFile matrixout = new TextFile(outmat, TextFile.W);
        String[] elems = matrix.readLineElems(TextFile.tab);
        List<String> list = Arrays.asList(elems);
        HashSet<String> meh = new HashSet<String>();
        meh.addAll(list);


        ArrayList<Integer> colsToKeep = new ArrayList<Integer>();

        colsToKeep.add(0);
        for (int i = 1; i < elems.length; i++) {
            if (keep.contains(elems[i])) {
                colsToKeep.add(i);
            }
        }
        if (colsToKeep.size() == 1) {
            System.out.println("Error: none of the columns match!");
            System.exit(-1);
        } else {
            System.out.println(colsToKeep.size() + " columns match.");
        }
        int[] keeparr = Primitives.toPrimitiveArr(colsToKeep);
        int lnctr = 0;
        while (elems != null) {
            String[] out = new String[keeparr.length];
            for (int c = 0; c < keeparr.length; c++) {
                out[c] = elems[keeparr[c]];
            }

            matrixout.writeln(Strings.concat(out, Strings.tab));
            lnctr++;
            if (lnctr % 10 == 0) {
                System.out.print("\r" + lnctr + " lines parsed");
            }
            elems = matrix.readLineElems(TextFile.tab);
        }

        matrix.close();
        matrixout.close();
    }
}
