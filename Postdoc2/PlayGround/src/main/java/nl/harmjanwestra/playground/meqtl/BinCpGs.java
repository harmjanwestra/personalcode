package nl.harmjanwestra.playground.meqtl;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class BinCpGs {
	
	public static void main(String[] args) {
		String matrix = "D:\\Work\\TMP\\oanejolien\\overlapratio_table_Bonder_ImputedGapped_unique_final_cpgs+feats+types.csv";
		String efile = "D:\\Work\\Projects\\2018-eQTMPredict\\bonder\\2015_09_02_cis_eQTMsFDR0.0-CpGLevel-filter.txtsplit.txt";
		
		int col = 1;
		int nrbins = 10;
		
		BinCpGs c = new BinCpGs();
		try {
			c.run(efile, matrix, col, nrbins);


			matrix = "D:\\Work\\TMP\\oanejolien\\overlapratio_table_Westra_ImputedGapped_unique_final_cpgs_Bonder_feats+types.csv";
			efile = "D:\\Work\\Projects\\2017-11-eQTLMeta\\eqtm\\eQTLsFDR0.0-SNPLevel-flipped.txt.gzsplit.txt";
			c.run(efile, matrix, col, nrbins);

			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public void run(String efile, String matrix, int col, int nrbins) throws IOException {
		
		int[][] ctrs = new int[2][nrbins+1];
		
		QTLTextFile f = new QTLTextFile(efile, QTLTextFile.R);
		ArrayList<EQTL> meqtl = f.readList();
		f.close();
		
		HashMap<String, Boolean> directions = new HashMap<String, Boolean>();
		for (EQTL e : meqtl) {
			String str = e.getRsName() + "_" + e.getProbe();
			if (e.getZscore() < 0) {
				directions.put(str, false);
			} else {
				directions.put(str, true);
			}
		}
		
		TextFile tf = new TextFile(matrix, TextFile.R);
		tf.readLine();
		String ln = tf.readLine();
		while (ln != null) {
			String[] elems = ln.split("\t");
			String eqtm = elems[0];
			boolean dir = directions.get(eqtm);
			Double d = Double.parseDouble(elems[col]);
			int bin = (int) Math.floor(d * nrbins);
			if (bin > nrbins) {
				bin = nrbins;
			}
			if (dir) {
				ctrs[0][bin]++;
			} else {
				ctrs[1][bin]++;
			}
			ln = tf.readLine();
			
		}
		tf.close();
		
		for (int b = 0; b < nrbins+1; b++) {
			double perc = (double) b / nrbins;
			System.out.println(b + "\t" + perc + "\t" + ctrs[0][b] + "\t" + ctrs[1][b]);
		}
		
	}
}
