package nl.harmjanwestra.playground.conv;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;

public class SNPMapToBed {
	
	public static void main(String[] args) {
		if (args.length < 2) {
			System.out.println("Usage: SNPMappings.txt.gz outbed");
		} else {
			SNPMapToBed b = new SNPMapToBed();
			try {
				b.run(args[0], args[1]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		
	}
	
	public void run(String in, String out) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfo = new TextFile(out, TextFile.W);
		
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		while (elems != null) {
			tfo.writeln("chr" + elems[0] + "\t" + elems[1] + "\t" + (Integer.parseInt(elems[1]) + 1) + "\t" + ctr);
			ctr++;
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfo.close();
		
	}
}
