package nl.harmjanwestra.playground.cis;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class LDValidator {
	
	
	public static void main(String[] args) {
		LDValidator v = new LDValidator();
//		String in1 = args[0];
//		String out1 = args[1];
//
//		try {
//			v.plot(in1, out1);
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		} catch (DocumentException e) {
//			e.printStackTrace();
//		}
//		System.exit(0);
		
		String snpannot = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\GiantSNPMappings_Filtered_on_Maf0.001_Indels.txt.gz";
		String refgenome = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-eqtlgen\\eur\\chrCHR.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz";
		String filter = null;
		
		snpannot = args[0];
		refgenome = args[1];
		String ldfile = args[2];
		String output = args[3];
		if (args.length > 4) {
			filter = args[4];
		}
		
		try {

//			String ldfile = "D:\\ld\\ENSG00000223972.txt.gz";
//			String output = "D:\\ld\\ENSG00000223972-compTo1Kg.txt";
//			refgenome = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-eqtlgen\\all\\chrCHR.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz";
			v.comparewithreference(ldfile, snpannot, refgenome, filter, output);
			
			ldfile = "D:\\ld\\ENSG00000227232.txt.gz";
			output = "D:\\ld\\ENSG00000227232-compTo1Kg.txt";
//			v.comparewithreference(ldfile, snpannot, refgenome, filter, output);
			
			String file1 = "D:\\ld\\ENSG00000223972.txt.gz";
			String file2 = "D:\\ld\\ENSG00000227232.txt.gz";
			String out = "D:\\ld\\ENSG00000227232-vs-ENSG00000223972.txt";

//			v.compareLDFiles(file1, file2, out);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void compareLDFiles(String file1, String file2, String outf) throws IOException {
		HashMap<String, String> map = new HashMap<>();
		
		TextFile tf = new TextFile(file1, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String ln = Strings.concat(elems, Strings.tab, 2, elems.length);
			map.put(elems[0] + "-" + elems[1], ln);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile tf2 = new TextFile(file2, TextFile.R);
		TextFile out = new TextFile(outf, TextFile.W);
		tf2.readLine();
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String query = elems[0] + "-" + elems[1];
			if (map.containsKey(query)) {
				String ln = Strings.concat(elems, Strings.tab) + "\t" + map.get(query);
				out.writeln(ln);
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		out.close();
	}
	
	
	public void comparewithreference(String ldfile, String snpannotfile, String refgenome, String filter, String output) throws IOException {
		
		// read possible snps
		HashSet<String> snps = new HashSet<String>();
		
		TextFile tf = new TextFile(ldfile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		int c = 0;
		while (elems != null) {
			
			String snp1 = elems[0];
			String snp2 = elems[1];
			
			snps.add(snp1);
			snps.add(snp2);
			
			c++;
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println(snps.size() + " snps to read, " + c + " lines in total");
		// load annotations
		
		HashMap<String, Feature> snpannot = loadAnnotation(snps, snpannotfile);
		
		// load reference
		
		VCFTabix tabix = null;
		tf.open();
		
		
		TextFile out = new TextFile(output, TextFile.W);
		out.writeln(tf.readLine() + "\tReferenceRSq");
		elems = tf.readLineElems(TextFile.tab);
		String prevsnp = null;
		Feature prevsnpObj = null;
		
		VCFVariant prevsnp1var = null;
		boolean[] samplefilter = null;
		DetermineLD ldc = new DetermineLD();
		int q = 0;
		HashMap<String, VCFVariant> snpToVariant = new HashMap<>();
		while (elems != null) {
			
			String snp1 = elems[0];
			String snp2 = elems[1];
			
			Feature snpf1 = snpannot.get(snp1);
			Feature snpf2 = snpannot.get(snp2);
			if (snpf1 != null && snpf2 != null) {
				if (prevsnp == null || !snp1.equals(prevsnp)) {
					
					if (tabix == null || !prevsnpObj.getChromosome().equals(snpf1.getChromosome())) {
						String vcffile = refgenome.replaceAll("CHR", "" + snpf1.getChromosome().getNumber());
						tabix = new VCFTabix(vcffile);
						if (filter != null) {
							samplefilter = tabix.getSampleFilter(filter);
						}
					}
					
					prevsnpObj = snpf1;
					prevsnp = snp1;
					VCFVariant var1 = null;
					if (snpToVariant.get(snp1) == null) {
						var1 = tabix.getVariant(snpf1, samplefilter);
						snpToVariant.put(snp1, var1);
					} else {
						var1 = snpToVariant.get(snp1);
					}
					prevsnp1var = var1;
				}
				VCFVariant snp2var = null;
				if (snpToVariant.get(snp2) == null) {
					snp2var = tabix.getVariant(snpf2, samplefilter);
					snpToVariant.put(snp2, snp2var);
					
				} else {
					snp2var = snpToVariant.get(snp2);
				}
				
				
				if (prevsnp1var != null && snp2var != null) {
					Pair<Double, Double> ld = ldc.getLD(prevsnp1var, snp2var);
					double rsq = ld.getRight();
					out.writeln(Strings.concat(elems, Strings.tab) + "\t" + rsq);
				}
				
			}
			elems = tf.readLineElems(TextFile.tab);
			q++;
			if (q % 10000 == 0) {
				System.out.print("\r" + q + " lines processed out of " + c + " : " + ((double) q / c));
			}
		}
		System.out.println();
		System.out.println("Done");
		tf.close();
		out.close();
		
	}
	
	private HashMap<String, Feature> loadAnnotation(HashSet<String> snps, String snpannotfile) throws IOException {
		HashMap<String, Feature> data = new HashMap<>();
		TextFile tf = new TextFile(snpannotfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		int c = 0;
		System.out.println("SNP annotation: " + snpannotfile);
		while (elems != null) {
			
			String snp = elems[2];
			if (snps.contains(snp)) {
				Chromosome chr = Chromosome.parseChr(elems[0]);
				Integer pos = Integer.parseInt(elems[1]);
				Feature f = new Feature(chr, pos, pos);
				data.put(snp, f);
			}
			
			c++;
			if (c % 100000 == 0) {
				System.out.print("\r" + c + " lines parsed");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println("");
		System.out.println("Done");
		
		System.out.println(data.size() + " annotations loaded");
		
		return data;
	}
	
	public void plot(String in, String out) throws IOException, DocumentException {
		
		Grid grid = new Grid(500, 500, 1, 1, 50, 50);
		
		ArrayList<Double> x = new ArrayList<Double>();
		ArrayList<Double> y = new ArrayList<Double>();
		
		TextFile tf = new TextFile(in, TextFile.R);
		String[] header = tf.readLineElems(TextFile.tab);
		int xpos = header.length - 2;
		int ypos = header.length - 1;
		System.out.println("X= " + header[xpos] + "\ty: " + header[ypos]);
		
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		while (elems != null) {
			
			
			Double xval = Double.parseDouble(elems[xpos]);
			Double yval = Double.parseDouble(elems[ypos]);
			
			x.add(xval);
			y.add(yval);
			
			ctr++;
			if (ctr % 10000 == 0) {
				System.out.print(ctr + " linss processed.\r");
			} else if (ctr > 10000000) {
				break;
			}
			elems = tf.readLineElems(TextFile.tab);
			
		}
		tf.close();
		System.out.println();
		System.out.println("Done.");
		
		DensityGraphPanel p = new DensityGraphPanel(1, 1);
		p.setPlotElems(true, false);
		
		p.setLabels("eQTLGen", "1KG");
		p.setDataRange(new Range(0, 0, 1, 1));
		p.setData(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
		p.setPlotElems(true, false);
		
		PearsonsCorrelation c = new PearsonsCorrelation();
		double corr = c.correlation(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
		p.setTitle("Pearson correlation: " + corr);
		System.out.println("Pearson correlation: " + corr);
		grid.addPanel(p);
		grid.draw(out);
		
	}
}
