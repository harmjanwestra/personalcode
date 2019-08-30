package nl.harmjanwestra.playground.cis.GWAS;

import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class CorrelateGenes {


	public static void main(String[] args) {

		String file = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\eQTLDump-sort-FDR-stripped-MS.txt.gz";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\eQTLDump-sort-FDR-stripped-MS-correlationsWithGenes.txt";
		CorrelateGenes c = new CorrelateGenes();
		try {
			c.correlate(file, output, 5e-8);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private class BetaPair {
		int chr;
		int pos;
		String rs;
		double beta1;
		double beta2;
		boolean beta1sig;
		boolean beta2sig;
	}

	public void correlate(String input, String output, double gwasthreshold) throws IOException {
		TextFile tf = new TextFile(input, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, ArrayList<BetaPair>> sets = new HashMap<>();
		int nrread = 0;
		while (elems != null) {

			BetaPair p = new BetaPair();
			p.beta1 = Double.parseDouble(elems[6]);
			p.beta2 = Double.parseDouble(elems[10]);
			p.rs = elems[2];
			p.chr = Integer.parseInt(elems[0]);
			p.pos = Integer.parseInt(elems[1]);

			double fdr = Double.parseDouble(elems[9]);
			double gwasp = Double.parseDouble(elems[12]);
			if (fdr < 0.05) {
				p.beta1sig = true;
			}
			if (gwasp < gwasthreshold) {
				p.beta2sig = true;
			}
			String gene = elems[5];
			ArrayList<BetaPair> pairs = sets.get(gene);
			if (pairs == null) {
				pairs = new ArrayList<>();
			}
			pairs.add(p);
			sets.put(gene, pairs);
			nrread++;
			if (nrread % 100000 == 0) {
				System.out.println(nrread + " read, " + sets.size() + " genes sofar.");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile out = new TextFile(output, TextFile.W);
		out.writeln("Gene\tNShared\tR");
		for (String gene : sets.keySet()) {
			ArrayList<BetaPair> pairs = sets.get(gene);

			if (pairs.size() > 30) {
				double[] x = new double[pairs.size()];
				double[] y = new double[pairs.size()];

				int ctr = 0;

				boolean genesig = false;
				boolean locussig = false;

				for (BetaPair p : pairs) {
					x[ctr] = p.beta1;
					y[ctr] = p.beta2;
					if (p.beta1sig) {
						genesig = true;
					}
					if (p.beta2sig) {
						locussig = true;
					}
					ctr++;
				}

				if (genesig && locussig) {
					SpearmansCorrelation c = new SpearmansCorrelation();
					double r = c.correlation(x, y);
					out.writeln(gene + "\t" + x.length + "\t" + r);
				}
			}
		}
		out.close();
	}
}
