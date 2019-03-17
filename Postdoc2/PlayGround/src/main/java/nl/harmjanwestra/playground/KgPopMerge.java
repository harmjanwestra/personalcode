package nl.harmjanwestra.playground;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class KgPopMerge {
	
	public static void main(String[] args) {
		
		try {
			TextFile tf = new TextFile("D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p3v5a-superpopulations.txt", TextFile.R);
			
			String[] elems = tf.readLineElems(TextFile.tab);
			
			HashMap<String, String> poptosuper = new HashMap<>();
			while (elems != null) {
				String pop = elems[0];
				String pop2 = elems[2];
				
				poptosuper.put(pop, pop2);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			TextFile tf2 = new TextFile("D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p3v5a-sampleinfo.txt", TextFile.R);
			TextFile tf3 = new TextFile("D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p3v5a-sampleinfo-wsuperpop.txt", TextFile.W);
			tf3.writeln("SuperPop\t" + tf2.readLine());
			elems = tf2.readLineElems(TextFile.tab);
			while (elems != null) {
				String pop = elems[1];
				String pop2 = poptosuper.get(pop);
				tf3.writeln(pop2 + "\t" + Strings.concat(elems, Strings.tab));
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();
			tf3.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
}
