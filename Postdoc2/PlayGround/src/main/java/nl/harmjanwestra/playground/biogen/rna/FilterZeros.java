package nl.harmjanwestra.playground.biogen.rna;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;

public class FilterZeros {

	public static void main(String[] args) {
		String input = "D:\\Sync\\SyncThing\\Data\\Ref\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR-CPM-TMM.txt.gz";
		String output = "D:\\Sync\\SyncThing\\Data\\Ref\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR-CPM-TMM-GenesWithMoreThan50PercZerosRemoved.txt.gz";
		FilterZeros z = new FilterZeros();
		try {
			z.run(input, output, 0.5);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void run(String in, String out, double percAllowed) throws IOException {
		TextFile outf = new TextFile(out, TextFile.W);
		TextFile tf = new TextFile(in, TextFile.R);
		outf.writeln(tf.readLine());

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			int nrz = 0;
			for (int i = 1; i < elems.length; i++) {
				double p = Double.parseDouble(elems[i]);
				if (p == 0.0d) {
					nrz++;
				}
			}
			double perc = (double) nrz / (elems.length - 1);
			if (perc < percAllowed) {
				outf.writeln(Strings.concat(elems, Strings.tab));
			}

			elems = tf.readLineElems(TextFile.tab);
		}
		outf.close();
		tf.close();
	}
}
