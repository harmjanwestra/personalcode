package nl.harmjanwestra.playground.conv;

import nl.harmjanwestra.utilities.enums.Chromosome;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class UpdateSNPMap {
	
	
	public static void main(String[] args) {
		UpdateSNPMap s = new UpdateSNPMap();
		
		try {
			s.runSNPMap(args[0], args[1]);
//			s.bimToBed(args[0], args[1], args[2]);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String snpmap, String indir) throws IOException {
		
		HashMap<String, String> posToRS = new HashMap<>();
		HashMap<String, String> rsToPos = new HashMap<>();
		
		System.out.println("Parsing: " + snpmap);
		TextFile tf = new TextFile(snpmap, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		while (elems != null) {
			Chromosome chr = Chromosome.parseChr(elems[0]);
			String pos = elems[1];
			String rsid = elems[2];
			posToRS.put(chr.getNumber() + ":" + pos, rsid);
			rsToPos.put(rsid, chr.getNumber() + ":" + pos);
			
			if (ctr % 10000 == 0) {
				System.out.print(ctr + " snps loaded.\r");
			}
			ctr++;
			elems = tf.readLineElems(TextFile.tab);
		}
		System.out.println();
		tf.close();
		
		System.out.println(posToRS.size() + " positions loaded");
		System.out.println(rsToPos.size() + " rsids loaded");
		
		for (int i = 1; i < 23; i++) {
			
			TextFile tfs = new TextFile(indir + "/chr" + i + "/SNPs.txt.gz", TextFile.R);
			System.out.println("Processing: " + tfs.getFileName());
			TextFile tfso = new TextFile(indir + "/chr" + i + "/SNPs-update.txt.gz", TextFile.W);
			TextFile tfsm = new TextFile(indir + "/chr" + i + "/SNPMappings-update.txt.gz", TextFile.W);
			int updated = 0;
			int removed = 0;
			int total = 0;
			String ln = tfs.readLine();
			while (ln != null) {
				if (ln.startsWith("rs")) {
					if (!rsToPos.containsKey(ln)) {
						tfso.writeln(ln + "-retired");
						tfsm.writeln("0\t0\t" + ln + "-retired");
						removed++;
					}
				} else {
					String[] snpelems = ln.split(":");
					if (snpelems.length == 2) {
						Chromosome chr = Chromosome.parseChr(snpelems[0]);
						String rs = posToRS.get(chr.getNumber() + ":" + snpelems[1]);
						if (rs == null) {
							tfso.writeln(ln + "-retired");
							tfsm.writeln("0\t0\t" + ln + "-retired");
							removed++;
						} else {
							tfso.writeln(rs);
							tfsm.writeln(chr.getNumber() + "\t" + snpelems[1] + "\t" + rs);
							updated++;
						}
					} else {
						tfso.writeln(ln + "-retired");
						tfsm.writeln("0\t0\t" + ln + "-retired");
						removed++;
					}
					
				}
				total++;
				ln = tfs.readLine();
			}
			tfso.close();
			tfsm.close();
			tfs.close();
			
			System.out.println(total + " written. " + removed + " removed. " + updated + " updated.");
			
		}
	}
	
	public void runSNPMap(String snpmap, String indir) throws IOException {
		
		HashMap<String, String> posToRS = new HashMap<>();
		HashMap<String, String> rsToPos = new HashMap<>();
		
		System.out.println("Parsing: " + snpmap);
		TextFile tf = new TextFile(snpmap, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		while (elems != null) {
			Chromosome chr = Chromosome.parseChr(elems[0]);
			String pos = elems[1];
			String rsid = elems[2];
			posToRS.put(chr.getNumber() + ":" + pos, rsid);
			rsToPos.put(rsid, chr.getNumber() + ":" + pos);
			
			if (ctr % 10000 == 0) {
				System.out.print(ctr + " snps loaded.\r");
			}
			ctr++;
			elems = tf.readLineElems(TextFile.tab);
		}
		System.out.println();
		tf.close();
		
		System.out.println(posToRS.size() + " positions loaded");
		System.out.println(rsToPos.size() + " rsids loaded");
		
		for (int i = 1; i < 23; i++) {
			
			TextFile tfs = new TextFile(indir + "/chr" + i + "/SNPMappings.txt.gz", TextFile.R);
			System.out.println("Processing: " + tfs.getFileName());
			TextFile tfso = new TextFile(indir + "/chr" + i + "/SNPs-update.txt.gz", TextFile.W);
			TextFile tfsm = new TextFile(indir + "/chr" + i + "/SNPMappings-update.txt.gz", TextFile.W);
			int updated = 0;
			int removed = 0;
			int total = 0;
			String ln = tfs.readLine();
			while (ln != null) {
				String[] lnelems = ln.split("\t");
				ln = lnelems[2];
				if (ln.startsWith("rs")) {
					if (!rsToPos.containsKey(ln)) {
						tfso.writeln(ln + "-retired");
						tfsm.writeln("0\t0\t" + ln + "-retired");
						removed++;
					}
				} else {
					String[] snpelems = ln.split(":");
					if (snpelems.length == 2) {
						Chromosome chr = Chromosome.parseChr(snpelems[0]);
						String rs = posToRS.get(chr.getNumber() + ":" + snpelems[1]);
						if (rs == null) {
							tfso.writeln(ln + "-retired");
							tfsm.writeln("0\t0\t" + ln + "-retired");
							removed++;
						} else {
							tfso.writeln(rs);
							tfsm.writeln(chr.getNumber() + "\t" + snpelems[1] + "\t" + rs);
							updated++;
						}
					} else {
						tfso.writeln(ln + "-retired");
						tfsm.writeln("0\t0\t" + ln + "-retired");
						removed++;
					}
					
				}
				total++;
				ln = tfs.readLine();
			}
			tfso.close();
			tfsm.close();
			tfs.close();
			
			System.out.println(total + " written. " + removed + " removed. " + updated + " updated.");
			
		}
	}
	
	public void convertRsToPos(String snpfile, String annotfile, String output) {
	
	}
	
	public void bimToBed(String snpfile, String bim, String bed) throws IOException {
		TextFile tf = new TextFile(snpfile, TextFile.R);
		ArrayList<String> snps = tf.readAsArrayList();
		tf.close();
		
		TextFile tf2 = new TextFile(bim, TextFile.R);
		String[] elems = tf2.readLineElems(TextFile.tab);
		HashMap<String, String> rsToPos = new HashMap<>();
		
		while (elems != null) {
			Chromosome chr = Chromosome.parseChr(elems[0]);
			String rs = elems[1];
			String pos = elems[3];
			rsToPos.put(rs, "chr" + chr.getNumber() + "\t" + pos + "\t" + pos);
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		
		int found = 0;
		int notfound = 0;
		TextFile bedout = new TextFile(bed, TextFile.W);
		for (String snp : snps) {
			String[] split = snp.split(":");
			if (split.length == 2) {
				Chromosome chr = Chromosome.parseChr(split[0]);
				bedout.writeln("chr" + chr.getNumber() + "\t" + split[1] + "\t" + split[1] + "\t" + snp);
			} else {
				String outstr = rsToPos.get(snp);
				if (outstr != null) {
					bedout.writeln(outstr + "\t" + snp);
					found++;
				} else {
					notfound++;
					bedout.writeln("chr0\t0\t0\t" + split);
				}
			}
		}
		bedout.writeln();
		System.out.println(found + " snps found in bim. " + notfound + " not found");
	}
	
}
