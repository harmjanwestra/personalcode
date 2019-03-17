package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashSet;

public class FilterGeneticSim {
	
	public static void main(String[] args) {
		
		String in = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\geneticsimilarity\\geneticsimilaritypairsAboveCorrelationThreshold0.25.txt.gz";
		FilterGeneticSim s = new FilterGeneticSim();
		try {
			s.run(in, null);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void run(String in, String out) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> dups = new HashSet<String>();
		while (elems != null) {
			String sample1 = elems[0];
			String sample2 = elems[1];
			Double sim = Double.parseDouble(elems[3]);
			if (sim >= 0.25 && !sample1.equals(sample2)) {
				dups.add(sample2);
			}
			
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		for (String s : dups) {
			System.out.println(s);
		}
		
		
	}
	
	
}
