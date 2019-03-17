package nl.harmjanwestra.playground.list;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class Compare {
	
	public static void main(String[] args) {
		String l1 = "D:\\temp\\allmethsamples.txt";
		String l2 = "D:\\temp\\allexpsamples.txt";
		String gtein = "D:\\temp\\allMTE.txt";
		String gteout = "D:\\temp\\allMTE-fixed.txt";
		Compare p = new Compare();
		try {
			p.run(l1, l2, gtein, gteout);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String l1, String l2, String gtein, String gteout) throws IOException {
		TextFile tf = new TextFile(l1, TextFile.R);
		HashSet<String> unique1 = new HashSet<String>();
		unique1.addAll(tf.readAsArrayList());
		tf.close();
		
		TextFile tf2 = new TextFile(l2, TextFile.R);
		HashSet<String> unique2 = new HashSet<String>();
		unique2.addAll(tf2.readAsArrayList());
		tf2.close();
		
		System.out.println(unique1.size());
		System.out.println(unique2.size());
		
		TextFile tfout = new TextFile(gteout, TextFile.W);
		TextFile tfin = new TextFile(gtein, TextFile.R);
		String[] elems = tfin.readLineElems(TextFile.tab);
		while (elems != null) {
			if (unique1.contains(elems[0]) && unique2.contains(elems[1])) {
				tfout.writeln(Strings.concat(elems, Strings.tab));
			}
			elems = tfin.readLineElems(TextFile.tab);
		}
		tfin.close();
		tfout.close();
		
	}
}
