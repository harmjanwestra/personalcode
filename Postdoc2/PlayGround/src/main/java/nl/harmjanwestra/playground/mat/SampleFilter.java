package nl.harmjanwestra.playground.mat;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashSet;

public class SampleFilter {
	
	public static void main(String[] args) {
		String select = "D:\\EMPBugFix\\sampleselect.txt";
		String mat = "D:\\Sync\\SyncThing\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR-CPM-TMM.txt.gz";
		String matout = "D:\\EMPBugFix\\expression.txt.gz";
		
		HashSet<String> samplesToInclude = new HashSet<String>();
		try {
			TextFile tfa = new TextFile(select, TextFile.R);
			String ln = tfa.readLine();
			while (ln != null) {
				samplesToInclude.add(ln);
				ln = tfa.readLine();
			}
			tfa.close();
			
			TextFile tf = new TextFile(mat, TextFile.R);
			TextFile outf = new TextFile(matout, TextFile.W);
			String[] elems = tf.readLineElems(TextFile.tab);
			boolean[] includecol = new boolean[elems.length];
			includecol[0] = true;
			String header = elems[0];
			for (int i = 1; i < elems.length; i++) {
				if (samplesToInclude.contains(elems[i])) {
					includecol[i] = true;
					header += "\t" + elems[i];
				}
			}
			outf.writeln(header);
			
			elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				String out = elems[0];
				for (int e = 1; e < elems.length; e++) {
					if (includecol[e]) {
						out += "\t" + elems[e];
					}
				}
				outf.writeln(out);
				elems = tf.readLineElems(TextFile.tab);
			}
			outf.close();
			tf.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
