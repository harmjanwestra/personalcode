package nl.harmjanwestra.playground.cis.perm;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class DetermineMaxAbsZScore {

    public static void main(String[] args) {


        DetermineMaxAbsZScore d = new DetermineMaxAbsZScore();
        if (args.length < 3) {
            System.out.println("Usage: genelist zmatloc outfile");
        } else {
            try {
                d.run(args[0], args[1], args[2]);
            } catch (IOException e) {
                e.printStackTrace();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    // generate table with max permuted z-score per gene
    public void run(String genelist, String zmatloc, String out) throws Exception {
        TextFile tf = new TextFile(genelist, TextFile.R);
        ArrayList<String> genes = tf.readAsArrayList();
        tf.close();

        System.out.println(genes.size() + " gene IDs loaded from: " + genelist);

        TextFile tfo = new TextFile(out, TextFile.W);
        // read first available gene
        int c = 0;
        for (String gene : genes) {
            String geneToLoad = zmatloc + gene + "-zmat.txt.gz";
            if (Gpio.exists(geneToLoad) && c == 0) {
                DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(geneToLoad);
                String datasets = Strings.concat(ds.getRowObjects(), Strings.tab);
                tfo.writeln("-\t" + datasets);
                c++;
            }
        }


        AtomicInteger ctr = new AtomicInteger();
        AtomicInteger finalCtr = ctr;
        ProgressBar progressBar = new ProgressBar(genes.size(), "Filtering...");
        IntStream.range(0, genes.size()).parallel().forEach(geneId -> {
            String gene = genes.get(geneId);
            String geneToLoad = zmatloc + gene + "-zmat.txt.gz";
            String samplesizetoload = zmatloc + gene + "-nmat.txt.gz";
            if (Gpio.exists(geneToLoad)) {
                // format: dataset x snps
                DoubleMatrixDataset<String, String> ds = null;
                DoubleMatrixDataset<String, String> dsn = null;
                try {
                    ds = DoubleMatrixDataset.loadDoubleData(geneToLoad);
                    dsn = DoubleMatrixDataset.loadDoubleData(samplesizetoload);

                    double[] maxZs = new double[ds.rows()];
                    for (int row = 0; row < ds.rows(); row++) {
                        double maxZ = Double.NaN;

                        for (int col = 0; col < ds.columns(); col++) {
                            double val = ds.getElementQuick(row, col);
                            if (!Double.isNaN(val)) {
                                val = Math.abs(val);
                                if (Double.isNaN(maxZ) || val > maxZ) {
                                    maxZ = val;
                                }
                            }

                        }
                        maxZs[row] = maxZ;
                    }



                    tfo.writelnsynced(gene + "\t" + Strings.concat(maxZs, Strings.tab));
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (Exception e) {
                    e.printStackTrace();
                }


                int ctrtmp = finalCtr.getAndIncrement();
                progressBar.set(ctrtmp);


            } else {
                System.out.println("Cannot find gene: " + geneToLoad);
            }
        });

        tfo.close();


    }
}
