package nl.harmjanwestra.playground.biogen.covariates;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;

public class PlinkGenomePihatFilter {


	public static void main(String[] args) {

		try {
			String filesfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\genomegzfiles.txt";
			TextFile tf = new TextFile(filesfile, TextFile.R);
			String ln = tf.readLine();
			PlinkGenomePihatFilter p = new PlinkGenomePihatFilter();
			double threshold = 0.125;
			while (ln != null) {
				String outf = ln.replaceAll(".gz", "-pihat" + threshold + ".txt");
				p.run(ln, outf, threshold);
				ln = tf.readLine();
			}
			tf.close();

		} catch (IOException e) {

		}

	}

	public void run(String infile, String outfile, double threshold) throws IOException {
		TextFile tf = new TextFile(infile, TextFile.R);
		TextFile tfo = new TextFile(outfile, TextFile.W);

		String ln = tf.readLine();

		int lnctr = 0;
		int beyondthr = 0;
		while (ln != null) {

			while (ln.startsWith(" ")) {
				ln = ln.substring(1, ln.length());
			}
			while (ln.contains("  ")) {
				ln = ln.replaceAll("  ", " ");
			}
			String[] elems = ln.split(" ");
			String pihatstr = elems[9];
			if (lnctr == 0) {
				tfo.writeln(Strings.concat(elems, Strings.tab));
			} else {
				double pihatd = Double.parseDouble(pihatstr);

				if (pihatd >= threshold) {
					tfo.writeln(Strings.concat(elems, Strings.tab));
					beyondthr++;
				}
			}
			lnctr++;
			ln = tf.readLine();
		}
		tfo.close();
		tf.close();

		System.out.println(beyondthr + " > " + threshold + " in " + outfile);
	}

}
