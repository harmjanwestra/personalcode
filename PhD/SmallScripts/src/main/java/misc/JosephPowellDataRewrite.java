/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.Set;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class JosephPowellDataRewrite {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            JosephPowellDataRewrite jpw = new JosephPowellDataRewrite();
            String in = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData/ExpressionData.txt.QuantileNormalized.Log2Transformed.txt.gz"; //"/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/RawExpressionData/ExpressionData.QuantileNormalized.Log2Transformed.txt.gz";////
            String probeToIlmn = "/Volumes/iSnackHD/Data/Projects/JosephPowell/probeToIlmnHT12v3.txt";
            String out = "/Volumes/iSnackHD/Data/Projects/JosephPowell/2013-11-15-Replication2/GroningenData/BloodHT12ExpressionData-QNorm-Log2Transformed.txt"; //"/Volumes/iSnackHD/Data/Projects/JosephPowell/2013-11-15-Replication2/EGCUTData/BloodHT12ExpressionData-QNorm-Log2Transformed.txt";// ;

            String samplesToInclude = null;//"/Volumes/iSnackHD/Data/Projects/JosephPowell/Data2/EGCUT/EGCUTExpSamples.txt";

            jpw.run(in, probeToIlmn, out, samplesToInclude);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void run(String in, String probeToIlmn, String out, String samplesToIncludeStr) throws IOException {

        Set<String> q = null;
        if (samplesToIncludeStr != null) {
            TextFile tfs = new TextFile(samplesToIncludeStr, TextFile.R);
            q = tfs.readAsSet(0, TextFile.tab);
            tfs.close();
        }
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(in, null, q);
        ds.transposeDataset();

        HashMap<String, String> probeToProbe = new HashMap<String, String>();
        TextFile tf = new TextFile(probeToIlmn, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String probe = elems[0];
            String ilmn = elems[1];
            probeToProbe.put(probe, ilmn);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        /*
         FID IID ILMN_xxxx ILMN_xxxx ILMN_xxxx ...
         f1  i1  x         x         x         ...
         f2  i2  x         x         x         ...
         f3  i3  x         x         x         ...

         */
        TextFile outfile = new TextFile(out, TextFile.W);
        String header = "FID\tIID";
        boolean[] includeCol = new boolean[ds.nrCols];
        for (int col = 0; col < ds.nrCols; col++) {
            String probe = ds.colObjects.get(col);
            String ilmn = probeToProbe.get(probe);
            if (ilmn != null) {

                header += "\t" + ilmn;
                includeCol[col] = true;
            } else {
                includeCol[col] = false;
            }
        }
        outfile.writeln(header);

        ProgressBar pb = new ProgressBar(ds.nrRows);
        for (int row = 0; row < ds.nrRows; row++) {
            StringBuilder rowOutput = new StringBuilder();
            rowOutput.append(0).append("\t").append(ds.rowObjects.get(row));
            for (int col = 0; col < ds.nrCols; col++) {
                if (includeCol[col]) {
                    rowOutput.append("\t").append(ds.rawData[row][col]);
                }
            }
            outfile.writeln(rowOutput.toString());
            pb.set(row);
        }
        pb.close();

        outfile.close();
    }
}
