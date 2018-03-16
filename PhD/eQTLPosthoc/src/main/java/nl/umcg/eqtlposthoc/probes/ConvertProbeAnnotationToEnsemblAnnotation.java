/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.probes;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ConvertProbeAnnotationToEnsemblAnnotation {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	try{
	    TextFile tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation.txt", TextFile.R);


	    TextFile tf2 = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-OnlyEnsemblAnnotation.txt", TextFile.W);

	    String[] elems = tf.readLineElems(TextFile.tab);
	    while(elems!=null){
		// Probe   Chr     Pos     Count   Ensembls

		//3       CAGCCATCTCTGCAGTTCTCTCAGTGCAGGCAGTTCTTCCTCTCAGGCTG      1       117446481-117446530     TTF2    ENSG00000116830 1       117404472       117447010
		if(elems[5].equals("-")){
		    elems[5] = "";
		}
		elems[5] = elems[5].replace(",", " ");
		if(!elems[3].equals("-") && !elems[3].equals("-1")){
		    elems[3] = elems[3].split("-")[0];
		}
		if(!elems[2].equals("-1")){
		    tf2.writeln(elems[0]+"\t"+elems[2]+"\t"+elems[3]+"\t"+0+"\t"+elems[5]);
		}

		elems = tf.readLineElems(TextFile.tab);
	    }

	    tf.close();
	    tf2.close();
	} catch (IOException e){
	    e.printStackTrace();
	}
    }
}
