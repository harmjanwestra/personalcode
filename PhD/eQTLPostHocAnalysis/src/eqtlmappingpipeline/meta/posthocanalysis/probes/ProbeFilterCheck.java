/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.probes;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ProbeFilterCheck {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here
	try {
	    TextFile tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt", TextFile.R);
	    ArrayList<String> probes = tf.readAsArrayList();
	    HashSet<String> probesSelected = new HashSet<String>();
	    probesSelected.addAll(probes);
	    tf.close();

	    TextFile tf2 = new TextFile("/Volumes/BackupDisk/MetaAnalysisFinal/trans/2012-05-02-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-CISTRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/TestedProbes.txt", TextFile.R);
	    String ln = tf2.readLine();
	    int nrpresent = 0;
	    int nrnotpresent = 0;
	    HashSet<String> probesTested = new HashSet<String>();
	    while (ln != null) {

		probesTested.add(ln);
		if (probesSelected.contains(ln)) {
		    nrpresent++;
		} else {
		    nrnotpresent++;
		}
		ln = tf2.readLine();
	    }

	    tf2.close();

	    System.out.println(nrpresent+"\t"+nrnotpresent);

	    for(int i=0; i<probes.size(); i++){
		if(!probesTested.contains(probes.get(i))){
		    System.out.println(probes.get(i)+"\tnot tested!");
		}
	    }
	} catch (IOException e) {
	}
    }
}
