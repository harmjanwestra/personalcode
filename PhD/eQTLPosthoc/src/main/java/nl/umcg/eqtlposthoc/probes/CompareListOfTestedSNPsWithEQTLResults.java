/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.probes;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class CompareListOfTestedSNPsWithEQTLResults {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here


	try {
	    TextFile q1 = new TextFile("/Volumes/BackupDisk/MetaAnalysisFinal/trans/2012-05-02-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-CISTRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/eQTLs.txt", TextFile.R);
	    HashSet<String> allsnps = new HashSet<String>();
	    String[] elems = q1.readLineElems(TextFile.tab);
	    while (elems != null) {
		if (elems.length > 2) {
		    allsnps.add(elems[1]);
		}

		elems = q1.readLineElems(TextFile.tab);
	    }
	    q1.close();

	    TextFile q2 = new TextFile("/Volumes/BackupDisk/MetaAnalysisFinal/trans/2012-05-02-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-CISTRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/TestedSNPs.txt", TextFile.R);
	    String[] snps = q2.readAsArray();
	    q2.close();

	    int nr = 0;
	    for(String snp: snps){
		if(!allsnps.contains(snp)){
		    System.out.println("SNP not found:\t"+snp);
		    nr++;
		}
	    }

	    System.out.println(nr);

	} catch (IOException e) {
	}
    }
}
