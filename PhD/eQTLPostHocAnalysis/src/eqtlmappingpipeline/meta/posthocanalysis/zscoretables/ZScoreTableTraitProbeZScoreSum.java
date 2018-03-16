/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.zscoretables;

import eqtlmappingpipeline.meta.posthocanalysis.SortableSNP;
import java.io.IOException;
import java.util.*;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class ZScoreTableTraitProbeZScoreSum {

    public void run() throws IOException {
	GWASCatalog catalog = new GWASCatalog();
	catalog.read("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog.txt");

	HashSet<GWASSNP> traitSNPs = new HashSet<GWASSNP>();


	HashSet<GWASTrait> selectedTraits = new HashSet<GWASTrait>();
	selectedTraits.addAll(Arrays.asList(catalog.getTraitsForCertainKey("Body mass index")));

	traitSNPs.addAll(Arrays.asList(catalog.getSNPsForTraitContainingKey("Body mass index")));



	TextFile probeTranslation = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-ProbeTranslationTable+H8HT12Conversion.log_reannotatedHG18_96PercIdentity.txt", TextFile.R);

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

	ArrayList<String> ldprunedSNPs = new ArrayList<String>();

	DetermineLD ldcalc = new DetermineLD();
	for (int s = 0; s < snpsWithAnnotation.size(); s++) {

	    String snp1 = snpsWithAnnotation.get(s).getname();
	    byte chr1 = snpsWithAnnotation.get(s).getchr();

	    ldprunedSNPs.add(snp1);
	    Integer snpid1 = genotypeData.getSnpToSNPId().get(snp1);
	    if (snpid1 == null) {
		System.out.println("SNP1 not in dataset! " + snp1);
	    } else {
		SNP snpobj1 = genotypeData.getSNPObject(snpid1);

		loader.loadGenotypes(snpobj1);
		loader.loadDosage(snpobj1);
		for (int s2 = s + 1; s2 < snpsWithAnnotation.size(); s2++) {
		    String snp2 = snpsWithAnnotation.get(s2).getname();
		    byte chr2 = snpsWithAnnotation.get(s2).getchr();
		    if (chr2 != chr1) {
			s = s2;
			break;
		    } else {
			// determine LD
			Integer snpid2 = genotypeData.getSnpToSNPId().get(snp2);
			if (snpid2 == null) {
			    System.out.println("SNP2 not in dataset! " + snp2);
			} else {
			    SNP snpobj2 = genotypeData.getSNPObject(snpid2);

			    loader.loadGenotypes(snpobj2);
			    loader.loadDosage(snpobj2);

			    double r2 = ldcalc.getRSquared(snpobj1, snpobj2, genotypeData, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);

			    snpobj2.clearGenotypes();
			    if (r2 < 0.5) {
				s = s2;
				break;
			    }
			}

		    }
		}

		snpobj1.clearGenotypes();
	    }

	}

	System.out.println(ldprunedSNPs.size() + " LD pruned SNPs remaining");
	HashSet<String> ldPrunedSNPsHash = new HashSet<String>();
	ldPrunedSNPsHash.addAll(ldprunedSNPs);

	TextFile zscoremat = new TextFile("/Volumes/BackupDisk/MetaAnalysis/cistrans/2012-02-14-HvH+SHIP+Rottedam+Inchianti+DILGOM-CisEffectsNotRegressedOut+EGCUT+Groningen-CisEffectsRegressedOut-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/metazscoretable.txt.gz", TextFile.R);

	elems = zscoremat.readLineElems(TextFile.tab);
	ArrayList<String> probeNames = new ArrayList<String>();
	for (int i = 3; i < elems.length; i++) {
	    probeNames.add(elems[i]);
	}

	double[] zsums = new double[elems.length - 3];
	int[] zsumsnpctr = new int[elems.length - 3];

	elems = zscoremat.readLineElems(TextFile.tab);
	while (elems != null) {

	    String snp = elems[0];

	    String coding = elems[1];
	    String[] bases = coding.split("/");
	    boolean atorcgsnp = false;

	    String assessed = elems[2];

	    GWASSNP snpObj = catalog.getSnpToObj().get(snp);

//	    if (traitSNPs.contains(snpObj)) {
	    if (ldPrunedSNPsHash.contains(snp)) {
		if (bases[0].equals(BaseAnnot.getComplement(bases[1]))) {
		    System.out.println("This is an A/T or C/G SNP!");
		    atorcgsnp = true;
		}

		GWASTrait[] traits = snpObj.getAssociatedTraitsArray();
		GWASTrait assoc = null;
		String riskAllele = null;
		for (GWASTrait t : traits) {
		    if (selectedTraits.contains(t)) {
			String traitSpecificRiskAllele = snpObj.getRiskAllele(t);
			if (traitSpecificRiskAllele != null) {
			    if (riskAllele == null) {
				assoc = t;
				riskAllele = traitSpecificRiskAllele;
				System.out.println("Same risk alleles for disease?\t" + snp + "\t" + coding + "\t" + assessed + "\t" + riskAllele + "\t" + traitSpecificRiskAllele);
			    } else {

				if (riskAllele.equals(traitSpecificRiskAllele)) {
				    System.out.println("Same risk alleles for disease?\t" + snp + "\t" + coding + "\t" + assessed + "\t" + riskAllele + "\t" + traitSpecificRiskAllele);
				} else {
				    System.out.println("!!!-- Different risk alleles for similar disease?\t" + snp + "\t" + coding + "\t" + assessed + "\t" + riskAllele + "\t" + assoc.getName() + "\t" + traitSpecificRiskAllele + "\t" + t.getName());
				}
			    }
			}
		    }
		}

		boolean flip = false;
		if (riskAllele != null) {
		    if (riskAllele.equals(assessed)) {
			// no flip
		    } else {
			if (riskAllele.equals(BaseAnnot.getComplement(assessed))) {
			    // no flip
			} else {
			    flip = true;
			}
		    }
		    Integer snpChr = snpToChr.get(snp);
		    Integer snpChrPos = snpToChrPos.get(snp);

		    for (int i = 3; i < elems.length; i++) {
			String probeName = probeNames.get(i - 3);
			Integer probeChr = probeToChr.get(probeName);
			Integer probeChrPos = probeToChr.get(probeName);
			try {
			    Double z = Double.parseDouble(elems[i]);

			    if (flip) {
				z = -z;
			    }

			    if (probeChr != null && snpChr != null) {

				if (snpChr.equals(probeChr)) {

				    if (probeChrPos != null && snpChrPos != null) {

					if (Math.abs(probeChrPos - snpChrPos) >= 5000000) {
					    zsums[i - 3] += Math.abs(z);
					    zsumsnpctr[i - 3]++;
					}
				    }

				} else {
				    zsums[i - 3] += Math.abs(z);
				    zsumsnpctr[i - 3]++;
				}

			    } else {
//				System.out.println(probeChr+"\t"+snpChr+"\t"+snp);
			    }



			} catch (Exception e) {
//			    if (probeName.equals("68561")) {
//
//				System.out.println(probeChr + "\t" + probeChrPos + "\t" + elems[i]);
//				e.printStackTrace();
//			    }

			}

		    }

		}

	    }

	    elems = zscoremat.readLineElems(TextFile.tab);
	}
	zscoremat.close();


	for (int i = 0; i < zsums.length; i++) {
	    if (zsumsnpctr[i] > 0) {
		System.out.println(i + "\t" + probeNames.get(i) + "\t" + probeToHUGO.get(probeNames.get(i)) + "\t" + zsums[i] + "\t" + zsumsnpctr[i] + "\t" + (zsums[i] / zsumsnpctr[i]));
	    } else {
		System.out.println(i + "\t" + probeNames.get(i) + "\t" + probeToHUGO.get(probeNames.get(i)) + "\t" + zsums[i] + "\t" + zsumsnpctr[i] + "\t" + 0);
	    }

	}

    }

    public static void main(String[] args) {
	try {
	    ZScoreTableTraitProbeZScoreSum s = new ZScoreTableTraitProbeZScoreSum();
	    s.run();
	} catch (Exception e) {
	    e.printStackTrace();
	}
    }
}
