/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.probes;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ProbeAnnotationCompare {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here

	try {

	    TextFile annotation1file = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-14-ProbeTranslationTable+H8HT12Conversion.txt", TextFile.R);

	    int nrProbes = 0;
	    HashSet<String> probesWithAnnotation = new HashSet<String>();
	    HashSet<String> probesInHT12v3 = new HashSet<String>();
	    String[] elems1 = annotation1file.readLineElems(TextFile.tab);
	    while (elems1 != null) {
		if (elems1.length > 4) {
		    String probe = elems1[0];
		    String chr = elems1[2];
		    if (!chr.equals("-1") && !chr.equals("-")) {
			probesWithAnnotation.add(probe);
		    }

		    if(!elems1[5].equals("-")){
			probesInHT12v3.add(probe);
		    }
		    nrProbes++;
		}
		elems1 = annotation1file.readLineElems(TextFile.tab);
	    }
	    annotation1file.close();

	    double ratio = (double) probesWithAnnotation.size() / nrProbes;
	    System.out.println(probesWithAnnotation.size() + " probes with an annotation out of " + nrProbes + "\t" + ratio);

	    TextFile annotation2file = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation.txt", TextFile.R);

	    String[] elems2 = annotation2file.readLineElems(TextFile.tab);
	    int nrWithAnnotatioInBoth = 0;
	    int nrWithOutAnnotationInBoth = 0;
	    int nrWithAnnotationInFirstButNotInLast = 0;
	    int nrWithAnnotationInLastButNotFirst = 0;
	    int nrprobes2 = 0;
	    int nrannotatedin2 = 0;

	    TextFile out = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt", TextFile.W);

	    int nrInHT12WithAnnotationInBoth = 0;
	    while (elems2 != null) {
		if (elems2.length > 4) {

		    String probe = elems2[0];
		    String chr = elems2[2];
		    if (!chr.equals("-1") && !chr.equals("-")) {
			nrannotatedin2++;
			if (probesWithAnnotation.contains(probe)) {
			    nrWithAnnotatioInBoth++;
			    if(probesInHT12v3.contains(probe)){
				out.writeln(probe);
				nrInHT12WithAnnotationInBoth++;
			    }
			} else {
			    nrWithAnnotationInLastButNotFirst++;
			    System.out.println("NotInFirst\t" + probe);
			}

		    } else {
			if (probesWithAnnotation.contains(probe)) {
			    nrWithAnnotationInFirstButNotInLast++;
			    System.out.println("NotInLast\t" + probe);
			} else {
			    nrWithOutAnnotationInBoth++;
			}
		    }
		    nrprobes2++;
		}
		elems2 = annotation2file.readLineElems(TextFile.tab);
	    }
	    annotation2file.close();

	    out.close();
	    System.out.println(nrWithAnnotatioInBoth + "\t" + nrWithAnnotationInLastButNotFirst + "\t" + nrWithAnnotationInFirstButNotInLast + "\t" + nrWithOutAnnotationInBoth);

	    System.out.println(nrprobes2);

	    System.out.println(nrannotatedin2);
	    System.out.println(nrInHT12WithAnnotationInBoth);
	} catch (IOException e) {
	    e.printStackTrace();
	}
    }
}
