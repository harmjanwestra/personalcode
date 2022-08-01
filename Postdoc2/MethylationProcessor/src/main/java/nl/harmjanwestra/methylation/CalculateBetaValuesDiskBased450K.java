package nl.harmjanwestra.methylation;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRandomAccessReader;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

public class CalculateBetaValuesDiskBased450K {


    public void run(String inU, String inM, String out) {


        System.out.println(inU);
        System.out.println(inM);
        System.out.println(out);

        try {
            DoubleMatrixDatasetRowIterable itU = new DoubleMatrixDatasetRowIterable(inU);
            DoubleMatrixDatasetRandomAccessReader<String, String> readerM = new DoubleMatrixDatasetRandomAccessReader<String, String>(inM);

            System.out.println(readerM.rows() + " x " + readerM.cols());

            ArrayList<String> uSamples = new ArrayList<String>(itU.getRows());

            Set<String> mCols = readerM.getHashCols().keySet();
            Set<String> uCols = itU.getCols();

            ArrayList<String> mColList = new ArrayList<String>(mCols);
            ArrayList<String> uColList = new ArrayList<String>(uCols);
            HashMap<String, Integer> mColMap = new HashMap<String, Integer>();
            for (int i = 0; i < mColList.size(); i++) {
                String probe = mColList.get(i);
                mColMap.put(probe, i);
            }

            int[] probeMap = new int[itU.getNrCols()];
            for (int i = 0; i < uColList.size(); i++) {
                Integer id = mColMap.get(uColList.get(i));
                if (id == null) {
                    System.out.println("Error: probe " + uColList.get(i) + " not found in " + inM);
                    System.exit(-1);
                } else {
                    probeMap[i] = mColMap.get(uColList.get(i));
                }
            }

            int rowctr = 0;
            double[] bvals = new double[itU.getNrCols()];
            DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(uColList, out);
            ProgressBar pb = new ProgressBar(itU.getNrRows(), "Calculating Beta-values.");
            for (double[] sampleU : itU) {
                String sample = uSamples.get(rowctr);

                Integer sampleMid = readerM.getHashRows().get(sample);
                // System.out.println(sample + "\tU: " + rowctr + "\tM: " + sampleMid);
                double[] sampleM = readerM.getRow(sampleMid);
                for (int d = 0; d < bvals.length; d++) {
                    double bval = (Math.max(sampleM[probeMap[d]], 0) + 100) / ((Math.max(sampleM[probeMap[d]], 0) + 100) + Math.max(sampleU[d], 0) + 100);
                    bvals[d] = bval;
                }
                writer.append(bvals, sample);
                rowctr++;
                pb.set(rowctr);
            }
            pb.close();
            itU.close();
            readerM.close();
            writer.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
