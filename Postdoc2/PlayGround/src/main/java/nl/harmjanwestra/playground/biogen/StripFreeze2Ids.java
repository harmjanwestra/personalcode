package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

public class StripFreeze2Ids {

	public static void main(String[] args) {
		String efile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-06-04-results\\cis\\2019-06-01-EUR-Cis-1mb-cortex-eQTLsFDR0.05-ProbeLevel.txt.gz";
		String efileout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-06-04-results\\cis\\2019-06-01-EUR-Cis-1mb-cortex-eQTLsFDR0.05-ProbeLevel-rsids.txt.gz";

		try {
			TextFile tf = new TextFile(efile, TextFile.R);
			TextFile tfo = new TextFile(efileout, TextFile.W);
			tfo.writeln(tf.readLine());
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				String id = elems[1];
				String[] idelems = id.split(":");
				elems[1] = idelems[2];
				tfo.writeln(Strings.concat(elems, Strings.tab));

				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			tfo.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
