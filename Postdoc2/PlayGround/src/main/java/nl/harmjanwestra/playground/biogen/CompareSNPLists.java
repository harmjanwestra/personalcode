package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class CompareSNPLists {

	public static void main(String[] args) {
//		String mayo = "D:\\TMP\\MAYO.snps.txt.gz";
//		String mayo = "D:\\TMP\\MAYO.snps.txt.gz";

		String sra = "D:\\TMP\\SRA.snps.txt.gz";
		String mayo = "D:\\TMP\\TargetALS.txt.gz";


		CompareSNPLists s = new CompareSNPLists();

		try {
			s.compare(mayo, sra, false);
			s.compare(mayo, sra, true);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void compare(String in1, String in2, boolean split) throws IOException {

		HashSet<String> set1 = loadset(in1, split);
		HashSet<String> set2 = loadset(in2, split);

		int ctr = 0;
		for (String s : set1) {
			if (set2.contains(s)) {
				ctr++;
			}
		}

		System.out.println(set1.size() + "\t" + set2.size() + "\t" + ctr);

	}

	private HashSet<String> loadset(String in1, boolean split) throws IOException {

		HashSet<String> output = new HashSet<String>();
		TextFile tf = new TextFile(in1, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {

			if (split) {
				String[] elems = ln.split(":");


				output.add(Strings.concat(elems, Strings.semicolon, 0, 3));
			} else {
				output.add(ln);
			}

			ln = tf.readLine();
		}


		tf.close();

		return output;
	}
}
