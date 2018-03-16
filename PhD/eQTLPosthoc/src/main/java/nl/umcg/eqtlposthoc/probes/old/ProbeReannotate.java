/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.probes.old;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ProbeReannotate {
    public void reannotateprobesInEQTLFile(String probeannotation, String eQTLFile) throws IOException {
	TextFile pa = new TextFile(probeannotation, TextFile.R);
	HashMap<String, String> probeToProbe = new HashMap<String, String>();

	String[] paelems = pa.readLineElemsReturnReference(TextFile.tab);
	paelems = pa.readLineElemsReturnReference(TextFile.tab);
	while (paelems != null) {
	    String id1 = paelems[4];
//	    String id2 = paelems[17];
//	    if (id1.equals("-")) {
//		id1 = id2;
//	    }
	    probeToProbe.put(paelems[0], id1);

	    paelems = pa.readLineElemsReturnReference(TextFile.tab);
	}
	pa.close();



	TextFile tf = new TextFile(eQTLFile, TextFile.R);
	TextFile tfReannotated = new TextFile(eQTLFile + "Reannotated.txt", TextFile.W);
	String[] elems = tf.readLineElemsReturnReference(tf.tab);

	String output = elems[0] + "";
	for (int i = 1; i < elems.length; i++) {
	    output += "\t" + elems[i];
	}
	output += "\n";
	tfReannotated.write(output);


	elems = tf.readLineElemsReturnReference(tf.tab);
	while (elems != null) {
	    String probe = elems[7];
//            Integer probeId  = Integer.parseInt(probe);
	    elems[19] = probeToProbe.get(probe);
//            elems[5]         = ""+probeChr[probeId];
//            elems[6]         = ""+probeChrPos[probeId];
//            elems[16]        = probeSymbol[probeId];
//
	    output = elems[0] + "";
	    for (int i = 1; i < elems.length; i++) {
		output += "\t" + elems[i];
	    }
	    output += "\n";
	    tfReannotated.write(output);

	    elems = tf.readLineElemsReturnReference(tf.tab);
	}
	tf.close();
	tfReannotated.close();
    }
}
