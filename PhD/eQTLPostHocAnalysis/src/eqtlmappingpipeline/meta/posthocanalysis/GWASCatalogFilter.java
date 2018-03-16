/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis;

import java.io.IOException;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASTrait;

/**
 *
 * @author harmjan
 */
public class GWASCatalogFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	try {
	    GWASCatalog c = new GWASCatalog();
	    c.read("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog.txt");
	    GWASTrait[] traits = c.getTraits();
	    for (GWASTrait t : traits) {
		System.out.println(t.getName());
	    }
	} catch (IOException e) {
	}
    }
}
