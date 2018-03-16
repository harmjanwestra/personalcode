/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class GOAnnotator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	try {
	    TextFile tf = new TextFile("/Volumes/BackupDisk/PathwayData/GO_BP/TermDescriptionsBP.txt", TextFile.R);
	    HashMap<String, String> map = new HashMap<String, String>();
	    String[] elems = tf.readLineElems(TextFile.tab);

	    while (elems != null) {
		map.put(elems[0], elems[1]);

		elems = tf.readLineElems(TextFile.tab);
	    }

	    tf.close();

	    TextFile tf2 = new TextFile("/Volumes/BackupDisk/MetaAnalysis/cistrans/2012-02-14-HvH+SHIP+Rottedam+Inchianti+DILGOM-CisEffectsNotRegressedOut+EGCUT+Groningen-CisEffectsRegressedOut-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/ZScoreCorrelations/GO_BP-PathWayAnnotation-FDR0.05-Cutoff0.001558.txt", TextFile.R);
	    TextFile tf3 = new TextFile("/Volumes/BackupDisk/MetaAnalysis/cistrans/2012-02-14-HvH+SHIP+Rottedam+Inchianti+DILGOM-CisEffectsNotRegressedOut+EGCUT+Groningen-CisEffectsRegressedOut-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/ZScoreCorrelations/GO_BP-PathWayAnnotation-FDR0.05-Cutoff0.001558.txt+Desc.txt", TextFile.W);
	    String line = tf2.readLine();
	    while (line != null) {
		elems = line.split("\t");
		String desc = map.get(elems[5]);

		tf3.writeln(line + "\t" + desc);

		line = tf2.readLine();
	    }
	    tf2.close();

	    tf3.close();
	} catch (IOException e) {
	    e.printStackTrace();
	}
    }
}
