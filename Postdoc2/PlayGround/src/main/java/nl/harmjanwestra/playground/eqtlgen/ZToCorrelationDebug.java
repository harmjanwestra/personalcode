package nl.harmjanwestra.playground.eqtlgen;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;

import java.io.IOException;
import java.util.HashMap;

public class ZToCorrelationDebug {

	public static void main(String[] args) {

		String ludeFile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-09-03-PBMCvsWholeBlood\\lude\\TransEQTLEffectSizeDifferenceBetweenBloodAndPBMCsHT12V4-20191115.txt";
		String hjfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-09-03-PBMCvsWholeBlood\\HJ\\eQTLsFDR.txt.gz";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-09-03-PBMCvsWholeBlood\\lude\\TransEQTLEffectSizeDifferenceBetweenBloodAndPBMCsHT12V4-20191115-vsHJ.txt";
		String mastertable = "D:\\TMP\\rtest\\2018-07-26-transEQTL-Replication-MasterTable.txt.gz";
		ZToCorrelationDebug z = new ZToCorrelationDebug();
		try {
			z.run(ludeFile, hjfile, mastertable, output);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void run(String file1, String file2, String tablefile, String outfile) throws IOException {

		HashMap<String, Double> tablersq = new HashMap<String, Double>();
		int tablecol = 13;
		TextFile tfq = new TextFile(tablefile, TextFile.R);
		tfq.readLine();
		tfq.readLine();
		String[] elems1 = tfq.readLineElems(TextFile.tab);
		while (elems1 != null) {
			String eqtl = elems1[3] + "-" + elems1[0];
			try {
				double d = Double.parseDouble(elems1[tablecol]);
				tablersq.put(eqtl, d);
			} catch (NumberFormatException e) {

			}
			elems1 = tfq.readLineElems(TextFile.tab);
		}
		tfq.close();
		HashMap<String, Double> eqtlslude = new HashMap<>();

		TextFile tf = new TextFile(file1, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String eqtl = elems[0] + "-" + elems[1];
			double r = Double.parseDouble(elems[5]);
			eqtlslude.put(eqtl, r);

			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		tf = new TextFile(file2, TextFile.R);
		TextFile out = new TextFile(outfile, TextFile.W);
		out.writeln("EQTL\tLude\tHJ\tTable");
		tf.readLine();
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String eqtl = elems[1] + "-" + elems[4];

			String[] narr = elems[13].split(";");
			int n = 0;
			for (int i = 0; i < narr.length; i++) {
				try {
					int ds = Integer.parseInt(narr[i]);
					n += ds;
				} catch (NumberFormatException e) {

				}
			}

			Double z = Double.parseDouble(elems[10]);

			double r = ZScores.zToR(z, n);

			Double ludeR = eqtlslude.get(eqtl);
			Double tableR = tablersq.get(eqtl);

			if (ludeR != null && tableR != null) {
				out.writeln(eqtl + "\t" + ludeR + "\t" + r + "\t" + tableR);
			} else if (ludeR == null) {
				System.out.println(eqtl + " not found by Lude.");
			} else if (tableR == null) {
				System.out.println(eqtl + " not found in table");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		out.close();
	}
}
