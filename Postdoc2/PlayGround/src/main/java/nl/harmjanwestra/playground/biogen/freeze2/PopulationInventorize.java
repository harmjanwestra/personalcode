package nl.harmjanwestra.playground.biogen.freeze2;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class PopulationInventorize {

	public static void main(String[] args) {


		String[] files = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\AMPAD-MAYO-V2-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\AMPAD-MSBB-V2-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\AMPAD-ROSMAP-V2-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Bipseq_1M-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Bipseq_2pt5M-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Bipseq_5M-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Bipseq_h650-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Braineac-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\BrainGVEX-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\CMC_HBCC_set1-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\CMC_HBCC_set2-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\CMC_HBCC_set3-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\CMC-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\GTEx-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\GVEX-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\GVEXUpdate-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Integrative-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\LIBD_1M-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\LIBD_5M-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\LIBD_h650-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\NABEC-H550-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\NABEC-H610-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-SampleAssignment.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\UCLA_ASD-SampleAssignment.txt"
		};
		PopulationInventorize p = new PopulationInventorize();
		for (String f : files) {
			try {
				p.run(f);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	public void run(String file) throws IOException {

		TextFile tf = new TextFile(file, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);


		HashMap<String, Integer> p = new HashMap<String, Integer>();
		p.put("AFR", 0);
		p.put("AMR", 1);
		p.put("EAS", 2);
		p.put("EUR", 3);
		p.put("SAS", 4);


		int[] cts = new int[5];
		int total = 0;
		while (elems != null) {
			if (elems.length > 1) {
				String sample = elems[0];
				String pop = elems[1];
				int id = p.get(pop);
				cts[id]++;
				total++;
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(file + "\t" + total + "\t" + Strings.concat(cts, Strings.tab));

	}
}
