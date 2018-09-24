package nl.harmjanwestra.playground.cis;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class RewriteAlleleFrequencyFile {
	
	private HashMap<String, String> snpmap;
	
	public static void main(String[] args) {
		RewriteAlleleFrequencyFile aq = new RewriteAlleleFrequencyFile();
		
		String allelefile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-07-18-AlleleFrequencies\\2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF.txt.gz";
		String snpannot = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\GiantSNPMappings_Filtered_on_Maf0.001_Indels.txt.gz";
		String referenceAlleleFile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\snpsAndAlleles-uniq.txt.gz";
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-07-18-AlleleFrequencies\\2018-07-18-AlleleFrequenciesFlippedToAssessed.txt.gz";
		
		
		try {
//			aq.gtfToProbeAnnotationFile(allelefile, snpannot, referenceAlleleFile, out);
			String combofolder = "D:\\snpgenecombos\\";
			String genelist = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\sortedGeneSNPCombos.txt.gz-genes.txt";
			String outfolder = "D:\\allelefreqspergene\\";
			aq.splitPerGene(out, genelist, combofolder, outfolder);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	public void splitPerGene(String alleleFrequencyFile, String genelist, String genedeffolder, String out) throws IOException {
		TextFile tf = new TextFile(genelist, TextFile.R);
		ArrayList<String> allgenes = tf.readAsArrayList();
		tf.close();
		
		HashMap<String, Double> allelefreq = new HashMap<>();
		TextFile atf = new TextFile(alleleFrequencyFile, TextFile.R);
		String header = atf.readLine();
		String[] elems = atf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[0];
			Double q = Double.parseDouble(elems[1]);
			allelefreq.put(snp, q);
			elems = atf.readLineElems(TextFile.tab);
		}
		atf.close();
		
		for (String g : allgenes) {
			String genedeffile = genedeffolder + g + ".txt.gz";
			System.out.println(out + g + ".txt.gz");
			tf = new TextFile(genedeffile, TextFile.R);
			TextFile atfo = new TextFile(out + g + ".txt.gz", TextFile.W);
			atfo.writeln(header);
			
			elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				String snp = elems[0];
				if (allelefreq.containsKey(snp)) {
					atfo.writeln(snp + "\t" + allelefreq.get(snp));
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			atfo.close();
			
		}
		
		
	}
	
	public void run(String allelefile, String snpannot, String referenceAlleleFile, String out) throws IOException {
		
		System.out.println("Reading reference allele map: " + referenceAlleleFile);
		HashMap<String, String> referenceAlleleMap = new HashMap<String, String>();
		{
			TextFile tf = new TextFile(referenceAlleleFile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			int ctr = 0;
			while (elems != null) {
				String snp = elems[0];
				String alleles = elems[1] + ";" + elems[2];
				referenceAlleleMap.put(snp, Strings.cache(alleles));
				ctr++;
//                if (ctr % 10000 == 0) {
//                    System.out.print(ctr + " lines parsed sofar\r");
//                }
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			System.out.println();
		}
		System.out.println("Reference alleles loaded for: " + referenceAlleleMap.size() + " snps.");
		
		loadSNPMap(snpannot, referenceAlleleMap);
		
		
		TextFile tfout = new TextFile(out, TextFile.W);
		TextFile tfin = new TextFile(allelefile, TextFile.R);
		tfin.readLine();
		String[] elems = tfin.readLineElems(TextFile.tab);
		tfout.writeln("RSid\tAssessedAlleleFrequency");
		int written = 0;
		int nrlines = 0;
		
		while (elems != null) {
			// 1-16383448      A       C       6       120     372     0.867469879518072
			String snp = elems[0];
			String rs = snpmap.get(snp);
			if (rs == null) {
				rs = snp;
			}
			
			String allele1 = elems[1];
			String allele2 = elems[2];
			
			
			String refAlleles = referenceAlleleMap.get(rs);
			
			if (refAlleles != null) {
				String[] refAlleleElems = refAlleles.split(";");
				
				String alleles = allele1 + "/" + allele2;
				String alleleassessed = allele2;
				
				Boolean flip = BaseAnnot.flipalleles(refAlleleElems[0], refAlleleElems[1], alleles, alleleassessed);
				
				Double freq = Double.parseDouble(elems[6]);
				if (flip != null) {
					if (flip) {
						freq = 1 - freq;
					}
					
					tfout.writeln(rs + "\t" + freq);
					written++;
				}
			} else {
				System.out.println(Strings.concat(elems, Strings.tab));
			}
			
			nrlines++;
			elems = tfin.readLineElems(TextFile.tab);
		}
		System.out.println(written + " out of " + nrlines);
		tfout.close();
		tfin.close();
		
		
	}
	
	private void loadSNPMap(String snpmapfile, HashMap<String, String> snplist) throws IOException {
		
		
		System.out.println("Loading snp map: " + snpmapfile);
		TextFile tf = new TextFile(snpmapfile, TextFile.R);
		
		snpmap = new HashMap<String, String>();
		
		String[] elems = tf.readLineElems(TextFile.tab);
		int lnctr = 0;
		while (elems != null) {
			if (elems.length >= 3) {
				String snp = elems[2];
				
				if (snplist == null || snplist.containsKey(snp)) {
					snp = new String(elems[2].getBytes(), "UTF-8");
					Chromosome chr = Chromosome.parseChr(elems[0]);
					Integer pos = Integer.parseInt(elems[1]);
					
					String id = chr.getNumber() + "-" + pos;
					snpmap.put(id, snp);
					
				}
			}
			lnctr++;
			if (lnctr % 10000 == 0) {
				System.out.print("\r" + lnctr + " lines parsed");
				
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println();
		
		System.out.println("SNP map has " + snpmap.size() + " elements.");
		
	}
}
