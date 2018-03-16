/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.qc;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ProbeMappingComparison {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here

	try {
	    TextFile tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/OldAnnotationfiles/2011-10-06-ProbeTranslationTable+H8HT12Conversion.log_reannotatedHG18_96PercIdentity.txt", TextFile.R);

	    String[] elems = tf.readLineElems(TextFile.tab);

	    HashSet<String> probesWithMapping = new HashSet<String>();
	    HashSet<String> probesWithOutMapping = new HashSet<String>();
	    while (elems != null) {

		if (!elems[2].equals("-")) {
		    probesWithMapping.add(elems[0]);
		} else {
		    probesWithOutMapping.add(elems[0]);
		}

		elems = tf.readLineElems(TextFile.tab);
	    }
	    tf.close();

	    System.out.println("W: "+ probesWithMapping.size());
	    System.out.println("W/O: "+ probesWithOutMapping.size());

	    tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-03-07-ProbeAnnotation96PercIdentity-ExonJunctionMappingsFixed-TranscriptMappingsFixed-ProbesWithWrongMappingLengthFilteredOut-SNPsPerProbe.txt", TextFile.R);

	    elems = tf.readLineElems(TextFile.tab);
	    int womapping = 0;
	    int wmapping = 0;
	    int overlap = 0;
	    int difference = 0;
	    HashSet<String> overlappingprobes = new HashSet<String>();
	    while (elems != null) {

		if (!elems[2].equals("-1")) {
		    wmapping++;
		    if (probesWithMapping.contains(elems[0])) {
			overlap++;
			overlappingprobes.add(elems[0]);
		    } else {

			difference++;
		    }
		} else {
		    womapping++;

		    if (probesWithOutMapping.contains(elems[0])) {
//			overlap++;
		    } else {
//			System.out.println("Probe has no mapping in new mapping file: "+elems[0]);
			difference++;
		    }
		}

		elems = tf.readLineElems(TextFile.tab);
	    }



	    tf.close();

	    TextFile tfout = new TextFile("", TextFile.W);
	    tfout.writeList(Arrays.asList(overlappingprobes.toArray(new String[0])));
	    tfout.close();

	    System.out.println("W: "+ wmapping);
	    System.out.println("W/O: "+ womapping);

	    System.out.println("");

	    System.out.println("Overlap: "+overlap);
	    System.out.println("Difference: "+difference);
	} catch (IOException e) {
	    e.printStackTrace();
	}
    }
}
