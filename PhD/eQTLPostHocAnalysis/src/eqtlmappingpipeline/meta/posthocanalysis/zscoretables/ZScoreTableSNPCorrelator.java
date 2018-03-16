/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.zscoretables;

import java.io.IOException;
import java.util.List;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class ZScoreTableSNPCorrelator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here
	// load the data...
	try {
	    TextFile output = new TextFile("/Data/PW/SNPCorrelations/RealData.txt", TextFile.W);
	    DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>("/Data/PW/metazscoretable.txt.gz-filtered.txt-ens.txt-collapsed.txt");


	    List<String> rows = ds.rowObjects;

	    String[] rownames = rows.toArray(new String[0]);

	    System.out.println("Starting correlation");
	    for(int i = 0; i<ds.nrRows; i++){
		System.out.println(i);
		double[] data = ds.rawData[i];
		for(int j=i+1; j<ds.nrRows; j++ ){
		    double[] data2 = ds.rawData[i];
		    double correlation = JSci.maths.ArrayMath.correlation(data, data2);
		    output.writeln(rownames[i]+"\t"+rownames[j]+"\t"+correlation+"\t"+(correlation*correlation));
		}
	    }

	    output.close();

	} catch (IOException e) {
	    e.printStackTrace();
	}
    }
}
