package nl.harmjanwestra.playground.cis;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;

public class eQTLRSquareCalc {
	
	
	public static void main(String[] args) {
		
		eQTLRSquareCalc r = new eQTLRSquareCalc();
		String in = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz";
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-RSq.txt.gz";
		
		try {
			r.run(in, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void run(String in, String out) throws IOException {
		
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfout = new TextFile(out, TextFile.W);
		String[] header = tf.readLineElems(TextFile.tab);
		tfout.writeln(Strings.concat(header, Strings.tab, 0, 12) + "\tR\tRSq");
		
		
		String[] elems = tf.readLineElems(TextFile.tab);
		
		int lnctr = 0;
		while (elems != null) {
			
			Integer n = Integer.parseInt(elems[12]);
			Double z = Double.parseDouble(elems[10]);
			
			double r = ZScores.zToR(z, n);
			double rsq = r * r;
			
			tfout.writeln(Strings.concat(elems, Strings.tab, 0, 12) + "\t" + r + "\t" + rsq);
			
			lnctr++;
			if (lnctr % 10000 == 0) {
				System.out.println(lnctr + " lines parsed.");
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfout.close();
		
	}
	
}

