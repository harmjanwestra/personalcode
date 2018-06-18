package nl.harmjanwestra.playground.conv;

import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

public class CleanExpTable {
	
	public static void main(String[] args) {
		CleanExpTable c = new CleanExpTable();
		try {
			c.run(args[0], args[1]);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void run(String in, String out) throws IOException {
		
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tf2 = new TextFile(out, TextFile.W);
		String[] elems = tf.readLineElems(Strings.tab);
		while (elems != null) {
			elems[1] = elems[1].split("\\.")[0];
			tf2.writeln(Strings.concat(elems, Strings.tab, 1, elems.length));
			elems = tf.readLineElems(Strings.tab);
		}
		tf2.close();
		tf.close();
	}
}
