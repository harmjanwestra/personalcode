/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ConvergenceCeliac {

    private static GWASCatalog catalog;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	try {
	    // TODO code application logic here
	    catalog = new GWASCatalog();

	    catalog.read("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog.txt");

	    GWASTrait[] traits = catalog.getTraits();
	    HashSet<String> traitSNPs = new HashSet<String>();
	    for (GWASTrait t : traits) {

		if (t.getName().toLowerCase().contains("celiac")) {
		    GWASSNP[] snps = t.getSNPs();
		    for (GWASSNP s : snps) {
			traitSNPs.add(s.getName());
			
			System.out.println(s.getName());
		    }
		}
	    }

	    TextFile tf = new TextFile("/Volumes/Data2/MetaAnalysisSoftware/ConvergentSNPsTrans0.5.txt", TextFile.R);

	    String[] elems = tf.readLineElems(TextFile.tab);
	    while (elems != null) {

		String gene = elems[3];
		String[] snpsassoc = elems[4].split(",");
		ArrayList<String> snpsselected = new ArrayList<String>();
		for (String s : snpsassoc) {
		    if (traitSNPs.contains(s)) {
			snpsselected.add(s);
		    }
		}

		if(snpsselected.size()>0){
		    for(String s: snpsselected){
			System.out.println(gene+"\t"+s);
		    }

		}


		elems = tf.readLineElems(TextFile.tab);
	    }

	} catch (IOException ex) {
	    Logger.getLogger(ConvergenceCeliac.class.getName()).log(Level.SEVERE, null, ex);
	}
    }
}
