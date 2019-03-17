package nl.harmjanwestra.playground.biogen.datasets;


import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class CMCLink {
	
	public static void main(String[] args) {
//        String gt = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\genotypepca\\CMC-Europeans.txt";
//        String rna = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\nonoutliers.txt";
//        String out = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\CMC-GTE-initial.txt";
		CMCLink c = new CMCLink();
		try {
//            c.run(gt, rna, out);
			
			String gte = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\CMC-GTE-initial.txt-filtered.txt";
			String rnaIDs = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\nonoutliers.txt";
			
			String out = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\CMC-GTE-initial.txt-filtered-rnaids.txt";
			
			c.replaceRNAIds(gte, rnaIDs, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void replaceRNAIds(String gte, String rnaIDs, String out) throws IOException {
		
		ArrayList<Pair<String, String>> pairs = new ArrayList<Pair<String, String>>();
		TextFile tf = new TextFile(gte, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			pairs.add(new Pair<String, String>(elems[0], elems[1]));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile tf2 = new TextFile(rnaIDs, TextFile.R);
		ArrayList<String> rnalist = tf2.readAsArrayList();
		HashSet<String> set = new HashSet<>();
		set.addAll(rnalist);
		tf2.close();
		
		TextFile tf3 = new TextFile(out, TextFile.W);
		for (Pair<String, String> p : pairs) {
			String g = p.getLeft();
			String r = p.getRight();
			
			boolean written = false;
			for (String s : rnalist) {
				
				if (s.contains(r) && !written) {
					tf3.writeln(g + "\t" + s);
					written = true;
				}
			}
			if (!written) {
				System.out.println("Could not find: " + r);
			}
			
		}
		tf3.close();
	}
	
	public void run(String gt, String rna, String out) throws IOException {
		TextFile tf = new TextFile(gt, TextFile.R);
		ArrayList<String> listgt = tf.readAsArrayList();
		HashSet<String> uniquegt = new HashSet<String>();
		uniquegt.addAll(listgt);
		tf.close();
		
		TextFile tf2 = new TextFile(rna, TextFile.R);
		ArrayList<String> listrna = tf2.readAsArrayList();
		HashSet<String> uniquerna = new HashSet<String>();
		uniquerna.addAll(listrna);
		tf2.close();
		
		HashSet<String> used = new HashSet<>();
		TextFile outf = new TextFile(out, TextFile.W);
		for (String g : uniquegt) {
			String[] split = g.split("_");
			String n = split[1];
			while (n.length() < 3) {
				n = "0" + n;
			}
			String combo = split[0] + "_" + n;
			
			if (g.equals("MSSM_3")) {
				System.out.println("?");
			}
			
			for (String s : uniquerna) {
				if (s.contains(combo)) {
					if (!used.contains(combo)) {
						outf.writeln(g + "\t" + s);
						used.add(s);
					} else {
						System.out.println("Dup: " + combo + "\t" + g);
					}
					
				}
			}
			if (!used.contains(combo)) {
				System.out.println("NF: " + combo + "\t" + g);
			}
		}
		outf.close();
	}
}
