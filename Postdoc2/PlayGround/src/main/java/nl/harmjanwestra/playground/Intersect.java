package nl.harmjanwestra.playground;


import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class Intersect {
	
	public static void main(String[] args) {
		String in1 = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\sortedGeneSNPCombos.txt.gz-genes.txt";
		String in2 = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-11-06-indepfx\\biosgeneswithindependentfx.txt";
		
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-11-06-indepfx\\biosgeneswithindependentfx-genesoverlapeqtlgen.txt";
		try {
			TextFile tf = new TextFile(in1, TextFile.R);
			ArrayList<String> list1 = tf.readAsArrayList();
			tf.close();
			
			HashSet<String> set1 = new HashSet<String>();
			set1.addAll(list1);
			
			
			TextFile tf2 = new TextFile(in2, TextFile.R);
			ArrayList<String> list2 = tf2.readAsArrayList();
			tf2.close();
			
			TextFile tf3 = new TextFile(out, TextFile.W);
			for (String s : list2) {
				if (set1.contains(s)) {
					tf3.writeln(s);
				}
			}
			tf3.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
