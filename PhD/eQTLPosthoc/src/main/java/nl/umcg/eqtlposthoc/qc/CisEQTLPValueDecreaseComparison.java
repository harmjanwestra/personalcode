/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.qc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class CisEQTLPValueDecreaseComparison {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here




	// take second best eQTL SNP for probe, considering their LD is ~0.9
	try {
	    // read SNP + probe combinations

	    HashMap<String, String> probeToSNPMap = new HashMap<String, String>();
	    HashMap<String, Double> probeToZScore = new HashMap<String, Double>();
	    TextFile eqtlfile = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-01-11-MetaAnalysis/cis/2012-02-22-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/QC/TruePositives.txteQTLProbesFDR0.05_sorted.txt", TextFile.R);

	    eqtlfile.readLine();

	    String[] elems = eqtlfile.readLineElems(TextFile.tab);
	    while (elems != null) {

		String snp = elems[eQTLTextFile.SNP];
		String probe = elems[eQTLTextFile.PROBE];

		probeToSNPMap.put(probe, snp);

		probeToZScore.put(probe, Math.abs(Double.parseDouble(elems[eQTLTextFile.METAZ])));
		elems = eqtlfile.readLineElems(TextFile.tab);
	    }
	    eqtlfile.close();



	    TriTyperGenotypeData ds = new TriTyperGenotypeData();
	    ds.load("/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/");

	    SNPLoader loader = ds.createSNPLoader();

	    DetermineLD ldcalc = new DetermineLD();

	    TextFile eqtlfile2 = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-01-11-MetaAnalysis/cis/2012-02-22-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/QC/TruePositives.txt_sorted.txt", TextFile.R);
	    eqtlfile2.readLine();
	    elems = eqtlfile2.readLineElems(TextFile.tab);
	    HashMap<String, String> probeToResult = new HashMap<String, String>();
	    int ctr = 0;
	    while (elems != null) {

		String snp = elems[eQTLTextFile.SNP];
		String probe = elems[eQTLTextFile.PROBE];

		String topSNP = probeToSNPMap.get(probe);
		if (topSNP != null && !probeToResult.containsKey(probe)) {

		    Double absZ = Math.abs(Double.parseDouble(elems[eQTLTextFile.METAZ]));
		    Double topZ = probeToZScore.get(probe);

		    Integer topsnpid = ds.getSnpToSNPId().get(topSNP);
		    if (topsnpid != null) {
			if (!snp.equals(topSNP) && absZ < topZ) {
			    Integer esnpid = ds.getSnpToSNPId().get(snp);
			    if (esnpid != null) {

				SNP topsnpobj = ds.getSNPObject(topsnpid);
				SNP esnpobj = ds.getSNPObject(esnpid);

				loader.loadGenotypes(esnpobj);
				loader.loadGenotypes(topsnpobj);

				double r2 = ldcalc.getRSquared(esnpobj, topsnpobj, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
				if (r2 > 0.80 && r2 < 0.85) {
				    // calculate the difference in Z-score.
				    double zDiff = topZ - absZ;
				    probeToResult.put(probe, snp + "\t" + r2 + "\t" + zDiff);
//				System.out.println(topSNP+"\t"+probe+"\t"+snp+ "\t"+r2+"\t"+zDiff);
				}

				topsnpobj.clearGenotypes();
				esnpobj.clearGenotypes();

			    }

			}
		    }
		}
		ctr++;
		if (ctr % 10000 == 0) {
		    System.out.println(ctr + " lines parsed");
		}
		elems = eqtlfile2.readLineElems(TextFile.tab);
	    }
	    eqtlfile2.close();

	    System.out.println("Done calculating");


	    TextFile tfout = new TextFile(eqtlfile.getFileName() + "-SNPsInLDWithTopSNP-ZScoreDifference-LD80-85.txt", TextFile.W);

	    // print results...
	    eqtlfile.open();
	    eqtlfile.readLine();
	    elems = eqtlfile.readLineElems(TextFile.tab);
	    while (elems != null) {

		String snp = elems[eQTLTextFile.SNP];
		String probe = elems[eQTLTextFile.PROBE];

		String result = probeToResult.get(probe);
//		System.out.println(snp+"\t"+probe+"\t"+result);

		tfout.writeln(snp + "\t" + probe + "\t" + result);

		elems = eqtlfile.readLineElems(TextFile.tab);
	    }
	    eqtlfile.close();
	    tfout.close();

	} catch (IOException e) {
	    e.printStackTrace();

	}
    }
}
