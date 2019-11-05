package nl.harmjanwestra.playground.biogen.eqtlposthoc;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

public class CountEQTSDiseases {


	public static void main(String[] args) {
		String file = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-23-PRS\\2019-10-23-PRS-eQTLsFDR0.05.txt";
		try {
			TextFile tf = new TextFile(file, TextFile.R);
			tf.readLine();

			HashMap<String, HashSet<String>> ctrs = new HashMap<String, HashSet<String>>();


			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {

				String[] diseaseElems = elems[1].split("\\_");
				String disease = Strings.concat(diseaseElems, Strings.dash, 0, diseaseElems.length - 1);
				String gene = elems[4];


				HashSet<String> ct = ctrs.get(disease);
				if (ct == null) {
					ct = new HashSet<String>();
				}
				ct.add(gene);
				ctrs.put(disease, ct);
				elems = tf.readLineElems(TextFile.tab);

			}
			tf.close();

			for (String key : ctrs.keySet()) {
				System.out.println(key + "\t" + ctrs.get(key).size());
			}

		} catch (IOException e) {
			e.printStackTrace();
		}


	}
}
