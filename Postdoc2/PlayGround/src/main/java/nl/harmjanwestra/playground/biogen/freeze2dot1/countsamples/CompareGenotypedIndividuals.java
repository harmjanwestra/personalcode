package nl.harmjanwestra.playground.biogen.freeze2dot1.countsamples;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class CompareGenotypedIndividuals {

	public static void main(String[] args) {
		String[] list = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\genotypeIndividuals\\BrainGVEX-V2-Individuals.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\genotypeIndividuals\\GVEXUpdate-Individuals.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\genotypeIndividuals\\GVEX-Individuals.txt"
		};

		for (int i = 0; i < list.length; i++) {
			try {
				TextFile tf = new TextFile(list[i], TextFile.R);
				HashSet<String> set1 = (HashSet<String>) tf.readAsSet(0, Strings.tab);
				tf.close();
				for (int j = i + 1; j < list.length; j++) {
					TextFile tf2 = new TextFile(list[j], TextFile.R);
					HashSet<String> set2 = (HashSet<String>) tf2.readAsSet(0, Strings.tab);
					tf2.close();
					int shared = 0;
					for (String s : set1) {
						if (set2.contains(s)) {
							shared++;
						}
					}
					System.out.println(i + "\t" + j + "\t" + set1.size() + "\t" + set2.size() + "\t" + shared);
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
}
