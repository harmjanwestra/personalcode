/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.qc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class DetermineTotalNumberOfeQTLProbesFromTwoFiles {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here

	try {
	    HashSet<String> uniqueProbes = new HashSet<String>();
	    TextFile tf = new TextFile("/Volumes/BackupDisk/MetaAnalysis/trans/2012-03-26-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-CISTRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/eQTLsFDR0.05.txt", TextFile.R);

	    String[] elems = tf.readLineElems(TextFile.tab);
	    while (elems != null) {
		uniqueProbes.add(elems[4]);
		elems = tf.readLineElems(TextFile.tab);
	    }
	    tf.close();

	    tf = new TextFile("/Volumes/BackupDisk/MetaAnalysis/trans/2011-12-21-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/eQTLsFDR0.05.txt", TextFile.R);
	    elems = tf.readLineElems(TextFile.tab);
	    while (elems != null) {
		uniqueProbes.add(elems[4]);
		elems = tf.readLineElems(TextFile.tab);
	    }
	    tf.close();


	    System.out.println(uniqueProbes.size());
	} catch (IOException e) {
	    e.printStackTrace();
	}
    }
}
