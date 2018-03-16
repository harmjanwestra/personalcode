/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.qc;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CisEQTLSampleCounter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here
	try{

	int[] datasetSizes = new int[9];
	datasetSizes[0] = 891;
	datasetSizes[1] = 963;
	datasetSizes[2] = 1240;
	datasetSizes[3] = 229;
	datasetSizes[4] = 762;
	datasetSizes[5] = 509;
	datasetSizes[6] = 611;
	datasetSizes[7] = 43;
	datasetSizes[8] = 63;


	TextFile tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-01-11-MetaAnalysis/cis/2012-02-22-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/QC/TruePositives.txteQTLProbesFDR0.05_sorted.txt", TextFile.R);

	String[] elems = tf.readLineElems(TextFile.tab);
	elems = tf.readLineElems(TextFile.tab);
	while(elems!=null){

	    String allzscores = elems[eQTLTextFile.DATASETZSCORE];
	    String[] zscores = Strings.semicolon.split(allzscores);

	    int numsamples = 0;
	    NumberFormatException e2 = null;
	    for(int i=0; i<zscores.length; i++){


		try{
		    Double.parseDouble(zscores[i]);
		    numsamples+=datasetSizes[i];
		} catch (NumberFormatException e){
		    e2 = e;
		}
	    }

	    System.out.println(numsamples);


	    elems = tf.readLineElems(TextFile.tab);
	}

	tf.close();
	} catch (IOException e){
	    e.printStackTrace();
	}
    }
}