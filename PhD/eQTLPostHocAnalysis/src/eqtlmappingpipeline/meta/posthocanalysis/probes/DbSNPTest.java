/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.probes;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.ncbi.dbsnp.SNPAnnotation;

/**
 *
 * @author harmjan
 */
public class DbSNPTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	DbSNPTest d = new DbSNPTest();
	try {
	    d.run();
	} catch (IOException ex) {
	    Logger.getLogger(DbSNPTest.class.getName()).log(Level.SEVERE, null, ex);
	}
    }

    private void run() throws IOException {

	SNPAnnotation s = new SNPAnnotation("/Data/SNPReferenceData/dbSNP/b130/b130_ContigInfo_36_3.bcp", "/Data/SNPReferenceData/dbSNP/b130/b130_SNPContigLoc_36_3.bcp", "reference");
	
    }


}
