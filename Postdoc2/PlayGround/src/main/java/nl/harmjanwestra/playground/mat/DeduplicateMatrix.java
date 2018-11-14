package nl.harmjanwestra.playground.mat;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class DeduplicateMatrix {
	
	
	public static void main(String[] args) {
		String in = "D:\\Sync\\SyncThing\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR.txt";
		String out = "D:\\Sync\\SyncThing\\geuvadis\\rnaseq-EUR\\GD660.GeneQuantCount-EUR-dedup.txt";
		
		DeduplicateMatrix d = new DeduplicateMatrix();
		try {
			d.run(in, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String in, String out) throws IOException {
		
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfout = new TextFile(out, TextFile.W);
		
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> visitedNames = new HashSet<String>();
		boolean[] inc = new boolean[elems.length];
		for (int i = 0; i < elems.length; i++) {
			if (!visitedNames.contains(elems[i])) {
				inc[i] = true;
				visitedNames.add(elems[i]);
			}
		}
		
		tfout.writeln(Strings.concat(elems, inc, Strings.tab));
		
		
		elems = tf.readLineElems(TextFile.tab);
		HashSet<String> visitedRows = new HashSet<String>();
		while (elems != null) {
			
			String id = elems[0];
			if (!visitedRows.contains(id)) {
				String outln = Strings.concat(elems, inc, Strings.tab);
				tfout.writeln(outln);
				visitedRows.add(id);
			}
			
			
			elems = tf.readLineElems(TextFile.tab);
		}
		
		tfout.close();
		tf.close();
		
	}
}
