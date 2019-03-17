package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class FilterEQTL {
	public static void main(String[] args) {
		
//		String efile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-21-lastminutegwas\\trans\\eQTLsFDR0.05.txt.gz";
//		String outdiseaseannot = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-21-lastminutegwas\\trans\\eQTLsFDR0.05-disease.txt";
		
		String efile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-21-lastminutegwas\\eQTLsFDR0.04.txt.gz";
		String outdiseaseannot = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-21-lastminutegwas\\eQTLsFDR0.04-diseaseannot.txt.gz";
		
		String[] disease = new String[]{
				"alzheimer",
				"ms",
				"parkinson"
		};
		String[] diseaseFile = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-21-lastminutegwas\\alzheimer.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-21-lastminutegwas\\ms.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-21-lastminutegwas\\parkinson.txt"
		};
		
		String outdir = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-21-lastminutegwas\\trans\\";
		String diseasefile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-21-lastminutegwas\\genetic_risk_factors_annotated_risk_allele_EBI_tested_traits_standardized_20180905.txt";
		
		FilterEQTL f = new FilterEQTL();
		for (int d = 0; d < disease.length; d++) {
			try {
//				f.run(diseaseFile[d], efile, outdir + disease[d] + ".txt");
				f.annotate(efile, diseasefile, outdiseaseannot);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		
	}
	
	private void annotate(String efile, String disease, String outdiseaseannot) throws IOException {
		
		HashMap<String, ArrayList<String>> diseasepersnp = new HashMap<String, ArrayList<String>>();
		
		TextFile tf = new TextFile(disease, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[0];
			String diseaseassoc = elems[8];
			
			ArrayList<String> list = diseasepersnp.get(snp);
			if (list == null) {
				list = new ArrayList<>();
			}
			list.add(diseaseassoc);
			diseasepersnp.put(snp, list);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile tfi = new TextFile(efile, TextFile.R);
		TextFile tfo = new TextFile(outdiseaseannot, TextFile.W);
		
		tfo.writeln(tfi.readLine() + "\tDisease");
		elems = tfi.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[1];
			ArrayList<String> list = diseasepersnp.get(snp);
			String diseasestr = "-";
			if (list != null) {
				diseasestr = Strings.concat(list, Strings.comma);
			}
			tfo.writeln(Strings.concat(elems, Strings.tab) + "\t" + diseasestr);
			elems = tfi.readLineElems(TextFile.tab);
		}
		
		tfo.close();
		tfi.close();
	}
	
	public void run(String snpfile, String efile, String out) throws IOException {
		TextFile tf = new TextFile(snpfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> snps = new HashSet<String>();
		while (elems != null) {
			String snp = elems[0];
			snps.add(snp);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile etf = new TextFile(efile, TextFile.R);
		String[] elems2 = etf.readLineElems(TextFile.tab);
		TextFile outf = new TextFile(out, TextFile.W);
		while (elems2 != null) {
			if (snps.contains(elems2[1])) {
				outf.writeln(Strings.concat(elems2, Strings.tab));
			}
			elems2 = etf.readLineElems(TextFile.tab);
		}
		etf.close();
		outf.close();
	}
}
