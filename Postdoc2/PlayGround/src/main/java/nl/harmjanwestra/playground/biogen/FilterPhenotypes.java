package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class FilterPhenotypes {

	public static void main(String[] args) {
		String f1 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-brain.phenotypes.txt";
		String f2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-RemainingSamplesAfterGenotypePCA.txt";
		String f3 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-brain.phenotypes-RemainingSamplesAfterGenotypePCA.txt";

		try {
			TextFile tf = new TextFile(f2, TextFile.R);
			HashSet<String> query = new HashSet<String>();
			query.addAll(tf.readAsArrayList());
			tf.close();
			TextFile tf2 = new TextFile(f1, TextFile.R);
			TextFile tf3 = new TextFile(f3, TextFile.W);
			tf3.writeln(tf2.readLine());
			String[] elems = tf2.readLineElems(TextFile.tab);
			while (elems != null) {
				if (query.contains(elems[1])) {
					tf3.writeln(Strings.concat(elems, Strings.tab));
				}
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();
			tf3.close();

		} catch (IOException e) {
			e.printStackTrace();
		}


	}
}
