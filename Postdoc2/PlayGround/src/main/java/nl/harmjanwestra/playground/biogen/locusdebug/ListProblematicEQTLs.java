package nl.harmjanwestra.playground.biogen.locusdebug;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class ListProblematicEQTLs {

	public static void main(String[] args) {
		String eqtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Cis-Genes\\Primary\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
		String faillist = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\2019-09-26-patch_genes-problematicsetcomparedtogtex.txt";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\eqtlfiltered\\";

		ListProblematicEQTLs e = new ListProblematicEQTLs();
		try {
			Gpio.createDir(output);
			e.run(faillist, eqtlfile, output);
		} catch (IOException ex) {
			ex.printStackTrace();
		}

	}

	public void run(String faillist, String eqtlfile, String output) throws IOException {

		HashSet<String> list = new HashSet<String>();


		TextFile tf = new TextFile(faillist, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			list.add(elems[0]);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile tf2 = new TextFile(eqtlfile, TextFile.R);
		String header = tf2.readLine();
		elems = tf2.readLineElems(TextFile.tab);
		TextFile outf1 = new TextFile(output + "ProblemEQTLs.txt", TextFile.W);
		TextFile outf2 = new TextFile(output + "NoProblemEQTLs.txt", TextFile.W);
		outf1.writeln(header);
		outf2.writeln(header);

		int ctr1 = 0;
		int ctr2 = 0;

		while (elems != null) {
			String gene = elems[4];
			if (list.contains(gene)) {
				outf1.writeln(Strings.concat(elems, Strings.tab));
				ctr1++;
			} else {
				outf2.writeln(Strings.concat(elems, Strings.tab));
				ctr2++;
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();

		System.out.println("Error: " + ctr1 + "\tNo error: " + ctr2);

		outf1.close();
		outf2.close();

	}
}
