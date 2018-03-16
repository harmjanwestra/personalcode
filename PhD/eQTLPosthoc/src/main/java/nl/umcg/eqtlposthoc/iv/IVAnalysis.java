/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.iv;

import java.io.IOException;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class IVAnalysis {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	try {

	    DoubleMatrixDataset<String,String> d = new DoubleMatrixDataset<String,String>("/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData4SNPPCRemoved/ExpressionData.txt.gz.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt");

	    Integer FADSId = (Integer) d.hashRows.get("2360020");

	    double[] dataForFADS1 = d.rawData[FADSId];

	    for (int i = 0; i < d.rowObjects.size(); i++) {

		if (!FADSId.equals(i)) {
		    double[] data2 = d.rawData[i];
		    double correlation = JSci.maths.ArrayMath.correlation(dataForFADS1, data2);
		    System.out.println(d.rowObjects.get(i) + "\t" + correlation + "\t" + (correlation * correlation));
		}
	    }

	} catch (IOException e) {
	    e.printStackTrace();
	}
    }
}
