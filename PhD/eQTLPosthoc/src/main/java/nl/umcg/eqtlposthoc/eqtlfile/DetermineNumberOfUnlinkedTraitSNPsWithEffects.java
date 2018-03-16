/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import nl.umcg.eqtlposthoc.SortableSNP;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

/**
 *
 * @author harmjan
 */
public class DetermineNumberOfUnlinkedTraitSNPsWithEffects {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here
	DetermineNumberOfUnlinkedTraitSNPsWithEffects d = new DetermineNumberOfUnlinkedTraitSNPsWithEffects();
	try {
	    d.run();
	} catch (IOException ex) {
	    Logger.getLogger(DetermineNumberOfUnlinkedTraitSNPsWithEffects.class.getName()).log(Level.SEVERE, null, ex);
	}
    }

    public void run() throws IOException {
	GWASCatalog catalog = new GWASCatalog();
	catalog.read("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog.txt");

	HashSet<GWASSNP> traitSNPs = new HashSet<GWASSNP>();


	HashSet<GWASTrait> selectedTraits = new HashSet<GWASTrait>();
	selectedTraits.addAll(Arrays.asList(catalog.getTraits()));

	traitSNPs.addAll(Arrays.asList(catalog.getSnps().toArray(new GWASSNP[0])));



	TextFile probeTranslation = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/OldAnnotationfiles/2011-10-06-ProbeTranslationTable+H8HT12Conversion.log_reannotatedHG18_96PercIdentity.txt", TextFile.R);

	HashMap<String, String> probeToHUGO = new HashMap<String, String>();
	HashMap<String, Integer> probeToChr = new HashMap<String, Integer>();
	HashMap<String, Integer> probeToChrPos = new HashMap<String, Integer>();

	String[] pbelems = probeTranslation.readLineElemsReturnObjects(TextFile.tab);
	while (pbelems != null) {

	    probeToHUGO.put(pbelems[0], pbelems[5]);
	    try {
		probeToChr.put(pbelems[0], Integer.parseInt(pbelems[2]));
		probeToChrPos.put(pbelems[0], Integer.parseInt(pbelems[3]));
	    } catch (NumberFormatException e) {
		System.out.println("Unable to parse chr pos: " + pbelems[0] + "\t" + pbelems[2] + "\t" + pbelems[4]);
	    }
	    pbelems = probeTranslation.readLineElemsReturnObjects(TextFile.tab);
	}

	probeTranslation.close();


	HashMap<String, Integer> snpToChr = new HashMap<String, Integer>();
	HashMap<String, Integer> snpToChrPos = new HashMap<String, Integer>();

	TextFile snpfile = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-27-SNPMappings-dbSNP130.txt.gz", TextFile.R);

	String[] elems = snpfile.readLineElems(TextFile.tab);
	while (elems != null) {
	    GWASSNP snpObj = catalog.getSnpToObj().get(elems[2]);

	    if (traitSNPs.contains(snpObj)) {
		System.out.println("Trait SNP detected: " + elems[2]);
		try {
		    snpToChr.put(elems[2], Integer.parseInt(elems[0]));
		    snpToChrPos.put(elems[2], Integer.parseInt(elems[1]));
		} catch (NumberFormatException e) {
		    System.out.println("Unable to parse chr pos: " + elems[0] + "\t" + elems[1] + "\t" + elems[2]);
		}
	    }
	    elems = snpfile.readLineElems(TextFile.tab);
	}
	snpfile.close();


	//
	GWASSNP[] traitSNPsArr = traitSNPs.toArray(new GWASSNP[traitSNPs.size()]);
	ArrayList<SortableSNP> snpsWithAnnotation = new ArrayList<SortableSNP>();
	for (int i = 0; i < traitSNPsArr.length; i++) {
	    String snpname = traitSNPsArr[i].getName();
	    Integer chr = snpToChr.get(snpname);
	    Integer chrPos = snpToChrPos.get(snpname);
	    if (chr != null && chrPos != null) {
		snpsWithAnnotation.add(new SortableSNP(snpname, chr.byteValue(), chrPos));
	    } else {
		System.out.println(snpname + " has null annotation: " + chr + "\t" + chrPos);
	    }
	}

	Collections.sort(snpsWithAnnotation);

	for (int i = 0; i < snpsWithAnnotation.size(); i++) {
	    System.out.println(snpsWithAnnotation.get(i).getname() + "\t" + snpsWithAnnotation.get(i).getchr() + "\t" + snpsWithAnnotation.get(i).getchrpos());
	}
//	System.exit(0);

	TriTyperGenotypeData genotypeData = new TriTyperGenotypeData();
	genotypeData.load("/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/");

	// prune the SNPs by LD

	SNPLoader loader = genotypeData.createSNPLoader();

//	ArrayList<String> ldprunedSNPs = new ArrayList<String>();
//
//	DetermineLD ldcalc = new DetermineLD();
//	for (int s = 0; s < snpsWithAnnotation.size(); s++) {
//
//	    String snp1 = snpsWithAnnotation.get(s).getname();
//	    byte chr1 = snpsWithAnnotation.get(s).getchr();
//
//	    ldprunedSNPs.add(snp1);
//	    Integer snpid1 = genotypeData.getSnpToSNPId().get(snp1);
//	    if (snpid1 == null) {
//		System.out.println("SNP1 not in dataset! " + snp1);
//	    } else {
//		SNP snpobj1 = genotypeData.getSNPObject(snpid1);
//
//		loader.loadGenotypes(snpobj1);
//		loader.loadDosage(snpobj1);
//		for (int s2 = s + 1; s2 < snpsWithAnnotation.size(); s2++) {
//		    String snp2 = snpsWithAnnotation.get(s2).getname();
//		    byte chr2 = snpsWithAnnotation.get(s2).getchr();
//		    if (chr2 != chr1) {
//			s = s2;
//			break;
//		    } else {
//			// determine LD
//			Integer snpid2 = genotypeData.getSnpToSNPId().get(snp2);
//			if (snpid2 == null) {
//			    System.out.println("SNP2 not in dataset! " + snp2);
//			} else {
//			    SNP snpobj2 = genotypeData.getSNPObject(snpid2);
//
//			    loader.loadGenotypes(snpobj2);
//			    loader.loadDosage(snpobj2);
//
//			    double r2 = ldcalc.getRSquared(snpobj1, snpobj2, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
//
//			    snpobj2.clearGenotypes();
//			    if (r2 < 0.5) {
//				s = s2;
//				break;
//			    }
//			}
//
//		    }
//		}
//
//		snpobj1.clearGenotypes();
//	    }
//
//	}
//
//	System.out.println(ldprunedSNPs.size() + " unlinked SNPs");

	HashMap<GWASSNP, HashSet<String>> cisEQTLsForSNPs = new HashMap<GWASSNP, HashSet<String>>();
	HashMap<GWASSNP, HashSet<String>> transEQTLsForSNPs = new HashMap<GWASSNP, HashSet<String>>();


	// old
	// "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-01-11-MetaAnalysis/trans/2011-12-21-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+HVH-CIS+TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/eQTLsFDR0.05.txt"
	// new
	//
	TextFile eqtlfile = new TextFile("/Volumes/BackupDisk/MetaAnalysisFinal/cistrans/2012-05-02-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut-TextOnlyOutput/eQTLsFDR0.05.txt", TextFile.R);

	HashSet<GWASSNP> selectedSNPsFromCatalog = new HashSet<GWASSNP>();
//	for (int i = 0; i < ldprunedSNPs.size(); i++) {
//	    selectedSNPsFromCatalog.add(catalog.getSnpToObj().get(ldprunedSNPs.get(i)));
//	}

	selectedSNPsFromCatalog.addAll(traitSNPs);

	elems = eqtlfile.readLineElems(TextFile.tab);
	int lnctr = 0;
	while (elems != null) {

	    String snp = elems[1];
	    String gene = elems[16];

	    GWASSNP s = catalog.getSnpToObj().get(snp);

//	    if (!gene.equals("-")) {
	    if (selectedSNPsFromCatalog.contains(s)) {

		boolean ciseffect = false;

		if (elems[2].equals(elems[5])) {
		    Integer snppos = Integer.parseInt(elems[3]);
		    Integer probepos = Integer.parseInt(elems[6]);

		    if (Math.abs(snppos - probepos) < 5000000) {
			ciseffect = true;
		    }

		}

		if (ciseffect) {
		    HashSet<String> cis = cisEQTLsForSNPs.get(s);
		    if (cis == null) {
			cis = new HashSet<String>();
		    }

		    cis.add(gene);

		    cisEQTLsForSNPs.put(s, cis);
		} else {
		    HashSet<String> trans = transEQTLsForSNPs.get(s);

		    if (trans == null) {
			trans = new HashSet<String>();
		    }

		    trans.add(gene);
		    transEQTLsForSNPs.put(s, trans);
		}




	    }
//	    }




	    elems = eqtlfile.readLineElems(TextFile.tab);
	    lnctr++;
	}

	System.out.println(lnctr + " eQTLs read in...");
	eqtlfile.close();

	GWASSNP[] snparr = selectedSNPsFromCatalog.toArray(new GWASSNP[selectedSNPsFromCatalog.size()]);
	int cistranseffects = 0;
	int transeffects = 0;
	int ciseffects = 0;

	int cisonly = 0;
	int transonly = 0;
	int total = 0;
	for (GWASSNP s : snparr) {

	    HashSet<String> cis = cisEQTLsForSNPs.get(s);
	    HashSet<String> trans = transEQTLsForSNPs.get(s);

	    if (cis != null) {
		ciseffects++;
	    }

	    if (trans != null) {
		transeffects++;
	    }

	    if (cis != null && trans != null) {
		cistranseffects++;
	    }

	    if (cis == null && trans != null) {
		transonly++;
	    }

	    if (cis != null && trans == null) {
		cisonly++;
	    }

	    if (cis != null || trans != null) {
		total++;
	    }
	}

	System.out.println("Number of snps with trans " + transeffects);
	System.out.println("Number of snps with cis " + ciseffects);
	System.out.println("Number of snps with cistrans " + cistranseffects);
	System.out.println("cis only " + cisonly);
	System.out.println("trans only " + transonly);
	System.out.println("total "+ total);
    }
}
