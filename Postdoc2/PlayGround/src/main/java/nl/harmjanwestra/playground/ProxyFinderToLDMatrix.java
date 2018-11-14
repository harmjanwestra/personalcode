package nl.harmjanwestra.playground;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import umcg.genetica.io.trityper.util.BaseAnnot;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class ProxyFinderToLDMatrix {
	
	public static void main(String[] args) {
		String reference = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-eqtlgen\\eur\\ldtest\\ids.txt";
		String in = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-eqtlgen\\eur\\ldtest\\ENS165092-Chr9_74605468-76605468.txt.gz";
		String out = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-eqtlgen\\eur\\ldtest\\ENS165092-Chr9_74605468-76605468-matrix-ldflippedTowardsEQTLGen.txt.gz";
		ProxyFinderToLDMatrix m = new ProxyFinderToLDMatrix();
		try {
			m.run(in, reference, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String list, String request, String out) throws IOException {
		
		HashSet<String> requestset = new HashSet<String>();
		HashMap<String, String> snpallelemap = new HashMap<String, String>();
		HashMap<String, String> snpalleleassessedmap = new HashMap<String, String>();
		
		TextFile tf1 = new TextFile(request, TextFile.R);
		String[] elems0 = tf1.readLineElems(TextFile.tab);
		while (elems0 != null) {
			String snp = elems0[0];
			String alleles = elems0[1];
			String alleleassessed = elems0[2];
			requestset.add(snp);
			snpallelemap.put(snp, alleles);
			snpalleleassessedmap.put(snp, alleleassessed);
			elems0 = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();
		System.out.println(requestset.size() + " SNPs requested ..");
		TextFile tf = new TextFile(list, TextFile.R);
		tf.readLine();
		HashMap<String, Integer> snpmap = new HashMap<String, Integer>();
		
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		ArrayList<String> snpsList = new ArrayList<String>();
		while (elems != null) {
			String snp1 = elems[2];
			String snp2 = elems[9];
			if (requestset.contains(snp1) && !snpmap.containsKey(snp1)) {
				snpmap.put(snp1, ctr);
				snpsList.add(snp1);
				ctr++;
			}
			if (requestset.contains(snp2) && !snpmap.containsKey(snp2)) {
				snpmap.put(snp2, ctr);
				snpsList.add(snp2);
				ctr++;
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(snpmap.size() + " snps found");
		
		double[][] matrix = new double[snpmap.size()][snpmap.size()];
		for (int d = 0; d < snpmap.size(); d++) {
			for (int d2 = d; d2 < snpmap.size(); d2++) {
				matrix[d][d2] = Double.NaN;
				matrix[d2][d] = Double.NaN;
			}
			matrix[d][d] = 1d;
		}
		
		tf.open();
		tf.readLine();
		elems = tf.readLineElems(TextFile.tab);
		
		int read = 0;
		int written = 0;
		while (elems != null) {
			String snp1 = elems[2];
			String alleles1 = elems[3].replaceAll(",", "/");
			String alleleAssessed1 = alleles1.split("/")[0];
			
			String snp2 = elems[9];
			String alleles2 = elems[10].replaceAll(",", "/");
			
			String alleleAssessed2 = alleles2.split("/")[0];
			
			Integer id1 = snpmap.get(snp1);
			Integer id2 = snpmap.get(snp2);
			
			if (id1 != null && id2 != null) {
				String allelesreq1 = snpallelemap.get(snp1);
				String allelesreq1ass = snpalleleassessedmap.get(snp1);
				String allelesreq2 = snpallelemap.get(snp2);
				String allelesreq2ass = snpalleleassessedmap.get(snp2);
				
				Boolean flip1 = BaseAnnot.flipalleles(allelesreq1, allelesreq1ass, alleles1, alleleAssessed1);
				Boolean flip2 = BaseAnnot.flipalleles(allelesreq2, allelesreq2ass, alleles2, alleleAssessed2);
				
				if (flip1 != null && flip2 != null) {
					boolean finalflip = false;
					if ((flip1 && flip2) || (!flip1 && flip2)) {
						finalflip = false;
					} else {
						finalflip = true;
					}
					
					double ld = Double.parseDouble(elems[elems.length - 2]);
					if (finalflip) {
						ld *= -1;
					}
					
					matrix[id1][id2] = ld;
					matrix[id2][id1] = ld;
					written++;
				}
				
			}
			read++;
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println(read + " read, " + written + " written.");
		
		String header = "-\t" + Strings.concat(snpsList, Strings.tab);
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln(header);
		for (int d = 0; d < matrix.length; d++) {
			String ln = snpsList.get(d) + "\t" + Strings.concat(matrix[d], Strings.tab);
			outf.writeln(ln);
		}
		outf.close();
		
	}
}
