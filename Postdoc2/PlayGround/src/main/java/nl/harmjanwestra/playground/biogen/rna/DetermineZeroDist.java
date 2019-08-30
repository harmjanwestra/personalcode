package nl.harmjanwestra.playground.biogen.rna;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

public class DetermineZeroDist {

	public static void main(String[] args) {
		String input = "D:\\Sync\\SyncThing\\Data\\Ref\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR-CPM-TMM.txt.gz";
		String output = "D:\\Sync\\SyncThing\\Data\\Ref\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR-CPM-TMM-nrZeros.txt";

		DetermineZeroDist z = new DetermineZeroDist();
		try {
			z.run(input, output);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void run(String in, String out) throws IOException {

		TextFile outf = new TextFile(out, TextFile.W);
		TextFile outf2 = new TextFile(out + "-Samples.txt", TextFile.W);
		TextFile tf = new TextFile(in, TextFile.R);
		String[] header = Strings.tab.split(tf.readLine());

		int[] nrZeroPerCol = new int[header.length - 1];
		int total = 0;
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			int nrz = 0;
			for (int i = 1; i < elems.length; i++) {
				double p = Double.parseDouble(elems[i]);
				if (p == 0.0d) {
					nrZeroPerCol[i - 1]++;
					nrz++;
				}
			}
			outf.writeln(elems[0] + "\t" + nrz + "\t" + (elems.length - 1));

			elems = tf.readLineElems(TextFile.tab);
			total++;
		}
		outf.close();

		for (int i = 1; i < header.length; i++) {
			outf2.writeln(header[i] + "\t" + nrZeroPerCol[i - 1] + "\t" + total);
		}

		outf2.close();
		tf.close();

	}

}
