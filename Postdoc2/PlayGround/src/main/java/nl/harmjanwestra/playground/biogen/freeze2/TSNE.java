package nl.harmjanwestra.playground.biogen.freeze2;

import com.jujutsu.tsne.TSneConfiguration;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import com.jujutsu.tsne.barneshut.BHTSne;
import com.jujutsu.tsne.barneshut.BarnesHutTSne;
import com.jujutsu.tsne.barneshut.ParallelBHTsne;
import com.jujutsu.utils.MatrixOps;
import com.jujutsu.utils.TSneUtils;

import java.util.ArrayList;
import java.util.HashSet;

public class TSNE {


    public static void main(String[] args) {
        String input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcacovariates\\EigenVectorsAfterCovarCorrection.txt.gz";
        String outplot = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcacovariates\\plots\\tsne.pdf";
        String outmat = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcacovariates\\plots\\tsne.txt.gz";
        TSNE t = new TSNE();
        try {
            t.run(input, outmat, outplot);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void run(String infile, String outmat, String outplot) throws Exception {


        HashSet<String> comps = new HashSet<String>();
        for (int i = 0; i < 50; i++) {
            comps.add("Comp" + i);
        }

        DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadSubsetOfTextDoubleData(infile, '\t', null, comps);


        int initial_dims = 55;
        double perplexity = 20.0;

        double[][] X = ds.getMatrix().toArray();
        X = MatrixOps.centerAndScale(X);
//        System.out.println(MatrixOps.doubleArrayToPrintString(X, ", ", 50, 10));
        BarnesHutTSne tsne;
        boolean parallel = true;
        if (parallel) {
            tsne = new ParallelBHTsne();
        } else {
            tsne = new BHTSne();
        }
        TSneConfiguration config = TSneUtils.buildConfig(X, 2, initial_dims, perplexity, 10000);
        double[][] Y = tsne.tsne(config);


        DoubleMatrixDataset dso = new DoubleMatrixDataset();
        dso.setMatrix(Y);
        ArrayList<String> rowids = new ArrayList<String>();
        ArrayList<String> colids = new ArrayList<String>();
        for (int i = 0; i < Y.length; i++) {
            rowids.add("Row" + i);
        }

        for (int j = 0; j < Y[0].length; j++) {
            colids.add("Col" + j);
        }
        dso.setColObjects(colids);
        dso.setRowObjects(rowids);
        dso.save(outmat);

        Grid grid = new Grid(300, 300, 1, 1, 100, 100);
        ScatterplotPanel p = new ScatterplotPanel(1, 1);
        p.setData(Y[0], Y[1]);
        p.setLabels("Dim1", "Dim2");
        grid.addPanel(p);
        grid.draw(outplot);
    }
}
