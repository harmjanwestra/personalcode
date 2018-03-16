/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.datasetcomparison;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.zip.DataFormatException;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultProbe;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;
import umcg.genetica.math.stats.Correlation;

/**
 *
 * @author harmjan
 */
public class CompareBinaryDatasets {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here

	try {

	    String[][] datasets = new String[9][2];
	    String[] platform = new String[9];
	    datasets[0][0] = "/Volumes/ADATA NH03/Marjolein/Results/2011-11-24-EGCUT-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemove/";
	    datasets[0][1] = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-10-03-EDCUT-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/";
	    platform[0] = "HT12v3";


	    datasets[1][0] = "/Volumes/ADATA NH03/Marjolein/Results/2011-11-09-SHIP_TREND-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/";
	    datasets[1][1] = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-19-03_SHIP-TREND-TRANS-40PCs-4GWASPCs-GeneticVextorsNotRemoved-CisEffectsRegressedOut/";
	    platform[1] = "HT12v3";


	    datasets[2][0] = "/Volumes/ADATA NH03/Marjolein/Results/2011-11-14-Groningen-BloodHT12-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/";
	    datasets[2][1] = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-03-09-Groningen-BloodHT12-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/";
	    platform[2] = "HT12v3";


	    datasets[3][0] = "/Volumes/ADATA NH03/Marjolein/Results/2011-11-14-Groningen-BloodH8v2-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/";
	    datasets[3][1] = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-03-09-Groningen-H8v2-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/";
	    platform[3] = "H8v2";


	    datasets[4][0] = "/Volumes/ADATA NH03/Marjolein/Results/2011-10-20-RotterdamStudy-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/";
	    datasets[4][1] = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-03-02-RotterdamStudy-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/";
	    platform[4] = "HT12v4";


	    datasets[5][0] = "/Volumes/ADATA NH03/Marjolein/Results/2011-11-15-DILGOM-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/";
	    datasets[5][1] = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-03-22-DILGOM-TRANS-40PCsRemoved-GeneticVectorsNotRemoved-4GWASPCsRemoved-CisEffectsRegressedOut/";
	    platform[5] = "HT12v3";


	    datasets[6][0] = "/Volumes/ADATA NH03/Marjolein/Results/2011-12-15-InChianti-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/";
	    datasets[6][1] = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-03-14-InChianti-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/";
	    platform[6] = "HT12v3";


	    datasets[7][0] = "/Volumes/ADATA NH03/Marjolein/Results/2011-12-10-HVH-HT12v3-TRANS-5PCs-4GWASPCs-GeneticVectorsNotRemoved/";
	    datasets[7][1] = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-03-19-HVH-v3-TRANS-5PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/";
	    platform[7] = "HT12v3";


	    datasets[8][0] = "/Volumes/ADATA NH03/Marjolein/Results/2011-12-22-HVH-HT12v4-TRANS-10PCs-4GWASPCs-GeneticVectorsNotRemoved/";
	    datasets[8][1] = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-03-19-HVH-v4-TRANS-10PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/";
	    platform[8] = "HT12v4";



	    HashSet<String> probesThatHaveBeenAlteredHT12v3 = CompareBinaryDatasets.getProbesFromFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-01-11-MetaAnalysis/cis/2012-02-22-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/eQTLProbesFDR0.05.txt-AnnotatedForHumanHT-12v3.txt");
	    HashSet<String> probesThatHaveBeenAlteredHT12v4 = CompareBinaryDatasets.getProbesFromFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-01-11-MetaAnalysis/cis/2012-02-22-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/eQTLProbesFDR0.05.txt-AnnotatedForHumanHT-12v4.txt");
	    HashSet<String> probesThatHaveBeenAlteredH8v2 = CompareBinaryDatasets.getProbesFromFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-01-11-MetaAnalysis/cis/2012-02-22-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/eQTLProbesFDR0.05.txtAnnotatedForH8v2.txt");



	    TextFile log = new TextFile("/Volumes/BackupDisk/comparisonlog.txt", TextFile.W);

	    for (int i = 0; i < datasets.length; i++) {

		log.writeln(ConsoleGUIElems.DOUBLELINE);
		log.writeln();
		log.writeln("Dataset1\t"+datasets[i][0]);
		log.writeln("Dataset2\t"+datasets[i][1]);

		BinaryResultDataset ds1 = new BinaryResultDataset(datasets[i][0], "Dataset", 0);

		BinaryResultDataset ds2 = new BinaryResultDataset(datasets[i][1], "Dataset", 0);


		BinaryResultProbe[] probes1 = ds1.getProbes();
		HashSet<String> probesInDs1 = new HashSet<String>();
		for (BinaryResultProbe p : probes1) {
		    probesInDs1.add(p.getName());
		}

		BinaryResultProbe[] probes2 = ds2.getProbes();
		int overlap = 0;
		for (BinaryResultProbe p2 : probes2) {
		    if (probesInDs1.contains(p2.getName())) {

			overlap++;
		    }
		}


		BinaryResultSNP[] snps1 = ds1.getSnps();
		HashSet<String> snpsInDs1 = new HashSet<String>();
		for (BinaryResultSNP s : snps1) {
		    snpsInDs1.add(s.getName());
		}

		BinaryResultSNP[] snps2 = ds2.getSnps();
		int overlapsnps = 0;
		for (BinaryResultSNP s2 : snps2) {
		    if (snpsInDs1.contains(s2.getName())) {
			overlapsnps++;
		    }
		}


		System.out.println("Probes in 1: " + probes1.length);
		System.out.println("Probes in 2: " + probes2.length);
		System.out.println("Overlap: " + overlap);

		log.writeln("Probes in 1: " + probes1.length);
		log.writeln("Probes in 2: " + probes2.length);
		log.writeln("Overlap: " + overlap);

		if (overlap < probes1.length || overlap < probes2.length) {
		    System.out.println("DIFFERENCE");
		}

		System.out.println("SNPs in 1: " + snps1.length);
		System.out.println("SNPs in 2: " + snps2.length);
		System.out.println("Overlap: " + overlapsnps);

		log.writeln("SNPs in 1: " + snps1.length);
		log.writeln("SNPs in 2: " + snps2.length);
		log.writeln("Overlap: " + overlapsnps);

		if (overlapsnps < snps1.length || overlapsnps < snps2.length) {

		    System.out.println("DIFFERENCE");

		}
		System.out.println("");

		log.writeln();

		ArrayList<String> probesToTest = new ArrayList<String>();
		HashSet<String> probesThatHaveBeenAltered = null;
		if (platform[i].equals("HT12v3")) {
		    probesThatHaveBeenAltered = probesThatHaveBeenAlteredHT12v3;
		} else if (platform[i].equals("HT12v4")) {
		    probesThatHaveBeenAltered = probesThatHaveBeenAlteredHT12v4;
		} else {
		    probesThatHaveBeenAltered = probesThatHaveBeenAlteredH8v2;
		}

		for (BinaryResultProbe p : probes1) {
		    if (!probesThatHaveBeenAltered.contains(p.getName())) {
			probesToTest.add(p.getName());
		    }
		}

		System.out.println("Dataset contains " + probesToTest.size() + " probes that have not had regression");
		log.writeln("Dataset contains " + probesToTest.size() + " probes that have not had regression");


		int probesThatShouldCorrelateButDont = 0;
		int probesThatShouldNotCorrelateButDo = 0;
		if (probesToTest.size() > 0) {
		    double[][] zscores1 = new double[probesToTest.size()][snps1.length];
		    for (int s = 0; s < snps1.length; s++) {
			Float[] zscores = ds1.readSNPZScores(snps1[s]);

			for (int p = 0; p < probesToTest.size(); p++) {
			    String name = probesToTest.get(p);
			    Integer id = ds1.getStringToProbe().get(name).getId();
			    zscores1[p][s] = zscores[id];
			}
		    }

		    int[] snpToSnp = new int[snps1.length];
		    for (int s = 0; s < snps1.length; s++) {
			BinaryResultSNP snp = ds2.getStringToSNP().get(snps1[s].getName());
			snpToSnp[s] = snp.getId();
//			System.out.println(s + "\t" + snpToSnp[s]);
		    }



		    double[][] zscores2 = new double[probesToTest.size()][snps1.length];
		    for (int s = 0; s < snps1.length; s++) {

			Float[] zscores = ds2.readSNPZScores(snps2[snpToSnp[s]]);
			for (int p = 0; p < probesToTest.size(); p++) {
			    String name = probesToTest.get(p);
			    Integer id = ds2.getStringToProbe().get(name).getId();
			    zscores2[p][s] = zscores[id];
			}
		    }

		    for (int p = 0; p < zscores2.length; p++) {

			double r = Correlation.correlate(zscores1[p], zscores2[p]);

			double r2 = r * r;
			if (r2 < 1) {
			    probesThatShouldCorrelateButDont++;
			    log.writeln(p+"\t"+probesToTest.get(p) + "\t" + r + "\t" + r2);
//			    System.out.println(probesThatHaveNotBeenAltered.get(p) + "\t" + r + "\t" + r2);
			}
//			System.out.println(probesThatHaveNotBeenAltered.get(p) + "\t" + r + "\t" + r2);
//			for(int s=0; s<snps1.length; s++){
//			    System.out.println(s+"\t"+zscores1[p][s]+"\t"+zscores2[p][s]);
//			}
		    }
		}



		probesToTest = new ArrayList<String>();
		for (BinaryResultProbe p : probes1) {
		    if (probesThatHaveBeenAltered.contains(p.getName())) {
			probesToTest.add(p.getName());
		    }
		}

		System.out.println("Dataset contains " + probesToTest.size() + " probes that have had regression");
		log.writeln("Dataset contains " + probesToTest.size() + " probes that have had regression");

		if (probesToTest.size() > 0) {
		    double[][] zscores1 = new double[probesToTest.size()][snps1.length];
		    for (int s = 0; s < snps1.length; s++) {
			Float[] zscores = ds1.readSNPZScores(snps1[s]);

			for (int p = 0; p < probesToTest.size(); p++) {
			    String name = probesToTest.get(p);
			    Integer id = ds1.getStringToProbe().get(name).getId();
			    zscores1[p][s] = zscores[id];
			}
		    }

		    int[] snpToSnp = new int[snps1.length];
		    for (int s = 0; s < snps1.length; s++) {
			BinaryResultSNP snp = ds2.getStringToSNP().get(snps1[s].getName());
			snpToSnp[s] = snp.getId();
//			System.out.println(s + "\t" + snpToSnp[s]);
		    }



		    double[][] zscores2 = new double[probesToTest.size()][snps1.length];
		    for (int s = 0; s < snps1.length; s++) {

			Float[] zscores = ds2.readSNPZScores(snps2[snpToSnp[s]]);
			for (int p = 0; p < probesToTest.size(); p++) {
			    String name = probesToTest.get(p);
			    Integer id = ds2.getStringToProbe().get(name).getId();
			    zscores2[p][s] = zscores[id];
			}
		    }

		    for (int p = 0; p < zscores2.length; p++) {

			double r = Correlation.correlate(zscores1[p], zscores2[p]);

			double r2 = r * r;
			if (r2 >= 1) {
			    probesThatShouldNotCorrelateButDo++;
			    log.writeln(p+"\t"+probesToTest.get(p) + "\t" + r + "\t" + r2);
//			    System.out.println(probesThatHaveNotBeenAltered.get(p) + "\t" + r + "\t" + r2);
			}
//			System.out.println(probesThatHaveNotBeenAltered.get(p) + "\t" + r + "\t" + r2);
//			for(int s=0; s<snps1.length; s++){
//			    System.out.println(s+"\t"+zscores1[p][s]+"\t"+zscores2[p][s]);
//			}
		    }
		}


		ds1.close();
		ds2.close();

		System.out.println("Probes that were regressed out, having a high correlation: " + probesThatShouldNotCorrelateButDo);
		System.out.println("Probes that were not regressed out, having a low correlation: " + probesThatShouldCorrelateButDont);

		System.out.println("");

		log.writeln("Probes that were regressed out, having a high correlation: " + probesThatShouldNotCorrelateButDo);
		log.writeln("Probes that were not regressed out, having a low correlation: " + probesThatShouldCorrelateButDont);

		log.writeln();
		System.out.println("");
	    }

	    log.close();
	} catch (IOException e) {

	    e.printStackTrace();
	} catch (DataFormatException e) {

	    e.printStackTrace();
	}
    }

    private static HashSet<String> getProbesFromFile(String string) throws IOException {
	TextFile tf = new TextFile(string, TextFile.R);

	HashSet<String> output = new HashSet<String>();
	String[] elems = tf.readLineElems(TextFile.tab);
	while (elems != null) {
	    output.add(elems[4]);
	    elems = tf.readLineElems(TextFile.tab);
	}
	tf.close();

	return output;
    }
}
