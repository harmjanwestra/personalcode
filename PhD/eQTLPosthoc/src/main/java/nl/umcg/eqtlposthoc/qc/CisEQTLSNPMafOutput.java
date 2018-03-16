/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.qc;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class CisEQTLSNPMafOutput {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here
//	try{
//
//	    TriTyperGenotypeData ds = new TriTyperGenotypeData();
//	    ds.load("/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/");
//
//	    SNPLoader loader = ds.createSNPLoader();
//
//	    TextFile tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-01-11-MetaAnalysis/cis/2012-02-22-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/QC/TruePositives.txteQTLProbesFDR0.05_sorted.txt", TextFile.R);
//
//	    tf.readLine();
//	    String[] elems = tf.readLineElems(TextFile.tab);
//	    while(elems != null){
//
//		String snp = elems[1];
//
//		Integer snpId = ds.getSnpToSNPId().get(snp);
//		if(snpId!=null){
//		    SNP snpobj = ds.getSNPObject(snpId);
//		    loader.loadGenotypes(snpobj);
//
//		    System.out.println(snpobj.getMAF());
//
//		    snpobj.clearGenotypes();
//		} else {
//		    System.out.println("NaN");
//		}
//
//
//		elems = tf.readLineElems(TextFile.tab);
//	    }
//
//	    tf.close();
//
//

	int ln = 0;
	try {

	    TextFile tf = new TextFile("/Volumes/ADATA NH03/Marjolein/Results/SNPMappings_SHIP-TREND.txt", TextFile.R);

	    HashMap<String, String> snpToR2 = new HashMap<String, String>();

	    String[] elems = tf.readLineElems(TextFile.tab);

	    while (elems != null) {

		String snp = elems[2];

		String r2Elems = elems[3];

		
		String[] pffrt = r2Elems.split(";");

		String[] qualstr = pffrt[1].split(":");

		String r2 = qualstr[1];

		snpToR2.put(snp, r2);

		elems = tf.readLineElems(TextFile.tab);
		ln++;
	    }
	    tf.close();


	    tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-01-11-MetaAnalysis/cis/2012-02-22-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/QC/TruePositives.txteQTLProbesFDR0.05_sorted.txt", TextFile.R);
//
	    tf.readLine();
	    elems = tf.readLineElems(TextFile.tab);
	    while (elems != null) {
		String snp = elems[1];
		String r2 = snpToR2.get(snp);
		System.out.println(r2);
		elems = tf.readLineElems(TextFile.tab);
	    }

	    tf.close();

	} catch (IOException e) {
	    e.printStackTrace();

	    System.out.println("at ln: "+ln);
	}


    }
}
