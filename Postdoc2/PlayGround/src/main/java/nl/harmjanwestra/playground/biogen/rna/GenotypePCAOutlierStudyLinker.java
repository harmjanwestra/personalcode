package nl.harmjanwestra.playground.biogen.rna;

import org.apache.poi.ss.formula.functions.T;
import org.w3c.dom.Text;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class GenotypePCAOutlierStudyLinker {

	public static void main(String[] args) {

		String allsamples = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\ENA\\allenasamples.txt";
		String samplelist = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\ENA\\Outliers.txt";
		String sampleToStudy = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\SampletoStudy.txt";

		try {
			HashMap<String, String> sampleStudyMap = new HashMap<String, String>();
			TextFile tf = new TextFile(sampleToStudy, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {

				sampleStudyMap.put(elems[0], elems[1]);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			HashSet<String> studies = new HashSet<String>();
			HashSet<String> visistedSamples = new HashSet<String>();
			TextFile tf2 = new TextFile(samplelist, TextFile.R);
			String ln = tf2.readLine();
			while (ln != null) {

				String study = sampleStudyMap.get(ln);
				System.out.println(ln + "\t" + study);
				studies.add(study);
				visistedSamples.add(ln);
				ln = tf2.readLine();
			}
			tf2.close();

			TextFile tf3 = new TextFile(allsamples, TextFile.R);
			ArrayList<String> allsamplelist = tf3.readAsArrayList();
			tf3.close();


			System.out.println();
			System.out.println("---");
			System.out.println("Studies");
			System.out.println("---");
			for (String s : studies) {
				System.out.println(s);
			}

			System.out.println();
			System.out.println("---");
			System.out.println("Missing samples");
			System.out.println("---");

			HashSet<String> studies2 = new HashSet<>();
			for (String s : allsamplelist) {
				if (!visistedSamples.contains(s)) {
					System.out.println(s + "\t" + sampleStudyMap.get(s));
					studies2.add(sampleStudyMap.get(s));
				}
			}
			System.out.println();
			System.out.println("---");
			System.out.println("Studies 2 ");
			System.out.println("---");
			for (String s : studies2) {
				System.out.println(s);
			}



		} catch (IOException e) {
			e.printStackTrace();
		}


	}
}
