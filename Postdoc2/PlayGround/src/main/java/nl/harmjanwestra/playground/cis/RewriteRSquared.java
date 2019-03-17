package nl.harmjanwestra.playground.cis;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;

public class RewriteRSquared {
	
	public static void main(String[] args) {
		RewriteRSquared r = new RewriteRSquared();
		try {
			r.run(args[0], args[1]);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String in, String out) throws IOException {
		
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfo = new TextFile(out, TextFile.W);
		tfo.writeln(tf.readLine());
		String[] elems = tf.readLineElems(TextFile.tab);
		
		while (elems != null) {
			Double r = Double.parseDouble(elems[elems.length - 1]);
			elems[elems.length - 2] = "" + r;
			elems[elems.length - 1] = "" + (r * r);
			tfo.writeln(Strings.concat(elems, Strings.tab));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfo.close();
		
	}
	
}

