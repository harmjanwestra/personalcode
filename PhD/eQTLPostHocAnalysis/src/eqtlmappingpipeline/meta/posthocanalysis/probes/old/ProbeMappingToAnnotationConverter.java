/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.probes.old;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ProbeMappingToAnnotationConverter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here

	try {



	    ProbeMappingToAnnotationConverter p = new ProbeMappingToAnnotationConverter();
	    p.readSequencesFromOldAnnotationFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/OldAnnotationfiles/2011-10-06-ProbeTranslationTable+H8HT12Conversion.log_reannotatedHG18_96PercIdentity.txt");
	    p.run("/Data/ProbeAnnotation/ProbeMappings/2012-02-21-IlluminaAll96PercentIdentity/all.intersect.paf_unique_combined", "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-02-22-ProbeAnnotation96PercIdentity-ExonJunctionMappingsFixed.txt");

	} catch (IOException e) {
	    e.printStackTrace();
	}
    }
    private HashMap<String, String> oldSeq;
    private ArrayList<String> probeNrs;

    public void run(String in, String out) throws IOException {

	TextFile tf = new TextFile(in, TextFile.R);

	TextFile tfout = new TextFile(out, TextFile.W);


	String[] elems = tf.readLineElems(TextFile.tab);

	HashMap<String, String> annotation = new HashMap<String, String>();

	while (elems != null) {

	    String probe = elems[1];
	    String chr = elems[3];
	    String chrstart = elems[4];

	    String chrend = elems[5];

	    String[] chrStartElems = chrstart.split(",");
	    String[] chrEndElems = chrend.split(",");


	    String finalAnnot = "";

	    for (int i = 0; i < chrEndElems.length; i++) {
		if (i == 0) {
		    finalAnnot += chrStartElems[i] + "-" + chrEndElems[i];
		} else {
		    finalAnnot += ":" + chrStartElems[i] + "-" + chrEndElems[i];
		}
	    }


	    String outStr = chr + "\t" + finalAnnot + "\t-";

	    annotation.put(probe, outStr);

	    elems = tf.readLineElems(TextFile.tab);
	}


	for (int i = 0; i < probeNrs.size(); i++) {
	    String probe = probeNrs.get(i);
	    String seq = oldSeq.get(probe);
	    String annot = annotation.get(probe);
	    if (annot == null) {
		annot = "-\t-\t-";
	    }
	    String outStr = probe + "\t" + seq + "\t" + annot;
	    tfout.writeln(outStr);
	}



	tfout.close();
	tf.close();

    }

    private void readSequencesFromOldAnnotationFile(String volumesData2MarjoleinHomeAccountmarjol) throws IOException {
	oldSeq = new HashMap<String, String>();

	probeNrs = new ArrayList<String>();
	TextFile old = new TextFile(volumesData2MarjoleinHomeAccountmarjol, TextFile.R);

	String[] elems = old.readLineElems(TextFile.tab);
	while (elems != null) {
	    oldSeq.put(elems[0], elems[1]);
	    probeNrs.add(elems[0]);
	    elems = old.readLineElems(TextFile.tab);
	}

	old.close();
    }
}
