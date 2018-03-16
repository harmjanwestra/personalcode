/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class AnnotateEQTLsWithSNPFunctionalPrediction {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here
	try {
	    TextFile eQTLFile = new TextFile("/Data/GeneticalGenomicsDatasets/HapMap-Stranger/Normalize/EQTLS/eQTLProbesFDR0.05.txt", TextFile.R);

	    HashSet<String> eSnps = new HashSet<String>();

	    String[] elems = eQTLFile.readLineElems(TextFile.tab);
	    elems = eQTLFile.readLineElems(TextFile.tab);
	    while (elems != null) {

		eSnps.add(elems[1]);
		elems = eQTLFile.readLineElems(TextFile.tab);
	    }

	    System.out.println(eSnps.size() + " eQTL SNPs loaded");
	    eQTLFile.close();

	    // now parse the snp reference data...

//	    // load the locus IDs corresponding to the reference sequence.
//	    TextFile dbsnpcontig = new TextFile("/Data/SNPReferenceData/dbSNP/b130/b130_ContigInfo_36_3_filtered.txt", TextFile.R);
//
//	    elems = dbsnpcontig.readLineElems(TextFile.tab);
//
//	    HashSet<String> contigs = new HashSet<String>();
//	    while (elems != null) {
//		contigs.add(elems[0]);
//		elems = dbsnpcontig.readLineElems(TextFile.tab);
//	    }
//	    dbsnpcontig.close();

	    // load snp annotaiton
	    TextFile dbsnp = new TextFile("/Data/SNPReferenceData/dbSNP/b130/b130_SNPContigLocusId_36_3.bcp", TextFile.R);

	    elems = dbsnp.readLineElems(TextFile.tab);
	    HashMap<String, Integer> functionalAnnotation = new HashMap<String, Integer>();

	    int ln = 0;
	    int nrannotated = 0;
	    while (elems != null) {


		if (eSnps.contains("rs" + elems[0])) {

		    // get functional class ID
		    String functionalClass = elems[11];
		    Integer numInClass = functionalAnnotation.get(functionalClass);
		    if (numInClass == null) {
			numInClass = 0;
		    }
		    numInClass++;
		    nrannotated++;

		    functionalAnnotation.put(functionalClass, numInClass);

		}
		elems = dbsnp.readLineElems(TextFile.tab);
		ln++;
		if (ln % 10000 == 0) {
		    System.out.println(ln + "\t" + nrannotated);
		}
	    }
	    dbsnp.close();

	    // load the annotation of the functional classes
	    TextFile functionalClassFile = new TextFile("/Data/SNPReferenceData/dbSNP/b130/SnpFunctionCode.bcp", TextFile.R);
	    elems = functionalClassFile.readLineElems(TextFile.tab);
	    while (elems != null) {
		String classId = elems[0];
		String classDesc = elems[1];
		Integer num = functionalAnnotation.get(classId);
		if (num == null) {
		    System.out.println("0\t(0%)\t" + classDesc);
		} else {
		    System.out.println(num + "/" + eSnps.size() + "\t(" + ((double) num / eSnps.size()) + "%)\t" + classDesc);
		}
		elems = functionalClassFile.readLineElems(TextFile.tab);
	    }

	    functionalClassFile.close();

	} catch (IOException e) {
	    e.printStackTrace();
	}
    }
}
