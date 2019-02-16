package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class QTLFilter {
	
	public static void main(String[] args) {
		
		String[] snpfiles = new String[]{
				"D:\\biogen\\GWASCatalog-2018-11-26-AD.txt",
				"D:\\biogen\\GWASCatalog-2018-11-26-ALS.txt",
				"D:\\biogen\\GWASCatalog-2018-11-26-MS.txt",
				"D:\\biogen\\GWASCatalog-2018-11-26-PD.txt"
		};
		
		QTLFilter f = new QTLFilter();
		
		String efile = "D:\\biogen\\10pcs-cis\\eQTLsFDR0.05-ProbeLevel.txt.gz";
		String[] out = new String[]{
				"D:\\biogen\\10pcs-cis\\eQTLsFDR0.05-ProbeLevel-AD.txt.gz",
				"D:\\biogen\\10pcs-cis\\eQTLsFDR0.05-ProbeLevel-ALS.txt.gz",
				"D:\\biogen\\10pcs-cis\\eQTLsFDR0.05-ProbeLevel-MS.txt.gz",
				"D:\\biogen\\10pcs-cis\\eQTLsFDR0.05-ProbeLevel-PD.txt.gz"
		};
		try {
			for (int e = 0; e < snpfiles.length; e++) {
				f.run(snpfiles[e], efile, out[e]);
				
			}
			
			efile = "D:\\biogen\\10pcs-trans\\eQTLsFDR0.05.txt.gz";
			out = new String[]{
					"D:\\biogen\\10pcs-trans\\eQTLsFDR0.05-AD.txt.gz",
					"D:\\biogen\\10pcs-trans\\eQTLsFDR0.05-ALS.txt.gz",
					"D:\\biogen\\10pcs-trans\\eQTLsFDR0.05-MS.txt.gz",
					"D:\\biogen\\10pcs-trans\\eQTLsFDR0.05-PD.txt.gz"
			};
			
			for (int e = 0; e < snpfiles.length; e++) {
				f.run(snpfiles[e], efile, out[e]);
				
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public void run(String snpfile, String efile, String out) throws IOException {
		HashSet<String> snps = new HashSet<String>();
		TextFile tf = new TextFile(snpfile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			snps.add(elems[2]);
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile outf = new TextFile(out, TextFile.W);
		TextFile inf = new TextFile(efile, TextFile.R);
		outf.writeln(inf.readLine());
		
		elems = inf.readLineElems(TextFile.tab);
		while (elems != null) {
			
			if (snps.contains(elems[1])) {
				outf.writeln(Strings.concat(elems, Strings.tab));
			}
			
			elems = inf.readLineElems(TextFile.tab);
		}
		inf.close();
		outf.close();
	}
	
	
}
