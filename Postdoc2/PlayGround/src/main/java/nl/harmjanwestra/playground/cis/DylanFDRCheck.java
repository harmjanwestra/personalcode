package nl.harmjanwestra.playground.cis;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.QTLTextFile;

import java.io.IOException;

public class DylanFDRCheck {
	
	
	public static void main(String[] args) {
		String meh = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz";
		
		DylanFDRCheck f = new DylanFDRCheck();
		try {
			f.run(meh);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void run(String efile) throws IOException {
		TextFile tf = new TextFile(efile, QTLTextFile.R);
		tf.readLine();
		String ln = tf.readLine();
		int ctr = 0;
		int ctr2 = 0;
		while (ln != null) {
			String[] elems = ln.split("\t");
			Double fdr = Double.parseDouble(elems[elems.length - 1]);
			if (fdr == 0) {
				ctr++;
			}
			if (fdr > 0) {
				break;
			}
			ln = tf.readLine();
		}
		tf.close();
		System.out.println(ctr + " eqtls!");
	}
}
