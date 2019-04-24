package nl.harmjanwestra.playground.biogen.rna;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class FilterRNAIds {

	public static void main(String[] args) {


		try {
			String link = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-brain.indtorna.txt";

			String outliers = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\ENA\\2019-04-12-GenotypePCAOutliers.txt";
			String remainingsamples = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-RemainingSamplesAfterGenotypePCA.txt";

			TextFile tf = new TextFile(outliers, TextFile.R);

			ArrayList<String> outlierlist = tf.readAsArrayList();

			tf.close();

			HashSet<String> outlierhash = new HashSet<String>();
			outlierhash.addAll(outlierlist);

			TextFile tfo = new TextFile(remainingsamples, TextFile.W);
			TextFile tf2 = new TextFile(link, TextFile.R);
			String[] elems = tf2.readLineElems(TextFile.tab);
			int ctr = 0;
			while (elems != null) {

				String ind = elems[0];
				String rna = elems[1];
				if (!outlierhash.contains(ind)) {
					tfo.writeln(rna);
				} else {
					ctr++;
				}

				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();
			tfo.close();
			System.out.println(ctr);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
