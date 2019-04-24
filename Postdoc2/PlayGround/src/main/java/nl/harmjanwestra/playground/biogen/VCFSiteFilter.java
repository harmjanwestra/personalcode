package nl.harmjanwestra.playground.biogen;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class VCFSiteFilter {

	public static void main(String[] args) {

	}

	public void run(String input, String sitelist, String chr, String output) throws IOException {

		Chromosome queryChr = Chromosome.parseChr(chr);
		HashSet<String> sites = new HashSet<>();
		TextFile tf = new TextFile(sitelist, TextFile.R);
		tf.readLine(); // skip header
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			Chromosome sitechr = Chromosome.parseChr(elems[0]);
			if (sitechr.equals(queryChr)) {
				String site = elems[0] + ":" + elems[1];
				sites.add(site);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(sites.size() + " sites from " + sitelist + " on chr " + chr);

		TextFile tfin = new TextFile(input, TextFile.R);
		TextFile tfout = new TextFile(input, TextFile.W);
		String ln = tfin.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {
				tfout.writeln(ln);
			} else {
				elems = Strings.subsplit(ln, Strings.tab,0,3);
				String query = elems[0] + ":" + elems[1];
				if (!sites.contains(query)) {
					tfout.writeln(ln);
				}
			}
			ln = tfin.readLine();
		}
		tfin.close();
		tfout.close();
	}
}
