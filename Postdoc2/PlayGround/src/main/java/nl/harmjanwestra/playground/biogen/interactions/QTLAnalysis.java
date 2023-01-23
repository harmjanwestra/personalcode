package nl.harmjanwestra.playground.biogen.interactions;

import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

public class QTLAnalysis {


    public static void main(String[] args) {


        String expdatafile = args[0]; //"/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/output-cortex/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.10PCAsOverSamplesRemoved.txt.gz";
        String gtfile = args[1];
        String gtefile = args[2];
        String genequery = args[3];
        String snpquery = args[4];


        try {
            HashSet<String> query = new HashSet<>();
            query.add(genequery);

            System.out.println("Loading expression data");
            DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadSubsetOfTextDoubleData(expdatafile, '\t', query, null);
            System.out.println(ds.rows() + " x " + ds.columns());
            HashMap<String, Double> gt = new HashMap<String, Double>();
            String alleles = null;
            TextFile tf = new TextFile(gtfile, TextFile.R);
            String[] header = tf.readLineElems(TextFile.tab);
            String[] elems = tf.readLineElems(TextFile.tab);
            System.out.println("Loading genotype data");
            while (elems != null) {
                if (elems[0].equals(snpquery)) {
                    alleles = elems[1];
                    for (int i = 3; i < elems.length; i++) {
                        gt.put(header[i], Double.parseDouble(elems[i]));
                    }
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile gtef = new TextFile(gtefile, TextFile.R);
            Map<String, String> gte = gtef.readAsHashMap(0, 1);
            gtef.close();

            ArrayList<Double> x = new ArrayList<>();
            ArrayList<Double> y = new ArrayList<>();

            for (String key : gte.keySet()) {
                String exp = gte.get(key);

                Double gtv = gt.get(key);
                Integer indid = ds.getColIndex(exp);
                Integer rowId = ds.getRowIndex(genequery);
                if (indid != null && rowId != null) {
                    Double expv = ds.getElement(rowId, indid);
                    if(gtv>-1) {
                        x.add(gtv);
                        y.add(expv);
                    }
                }
            }

            SpearmansCorrelation c = new SpearmansCorrelation();
            double r = c.correlation(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
            Correlation.correlationToZScore(x.size());
            double z = Correlation.convertCorrelationToZScore(x.size(), r);
            double p = ZScores.zToP(z);

            // correlate
            System.out.println(snpquery + "\t" + genequery + "\t" + alleles + "\t" + x.size() + "\t" + r + "\t" + z + "\t" + p);
        } catch (Exception e) {
            e.printStackTrace();
        }


    }
}
