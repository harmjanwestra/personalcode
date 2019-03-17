package nl.harmjanwestra.playground.biogen.datasets.gtex;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class MergePsychencode {
	
	public static void main(String[] args) {
		
		String efile = "D:\\eQTLTest\\DER-08b_hg19_eQTL.bonferroni.txt";
		String efileo = "D:\\eQTLTest\\DER-08b_hg19_eQTL.bonferroni-withalleles.txt";
		String allelefile = "D:\\eQTLTest\\SNP_Information_Table_with_Alleles.txt";
		
		
		MergePsychencode p = new MergePsychencode();
		try {
//			p.run(allelefile, efile, efileo);
			String file = "D:\\eQTLTest\\psychencode-DER-08a_hg19_eQTL.significant.txt";
			String out = "D:\\eQTLTest\\psychencode-DER-08a_hg19_eQTL.significant-topfx.txt";
			p.filtertopfx(file, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void filtertopfx(String file, String out) throws IOException {
		TextFile tf = new TextFile(file, TextFile.R);
		TextFile tfo = new TextFile(out, TextFile.W);
		
		String[] header = tf.readLineElems(TextFile.tab);
		int topsnpcol = 0;
		for (int i = 0; i < header.length; i++) {
			String v = header[i];
			if (v.equals("top_SNP")) {
				topsnpcol = i;
			}
		}
		
		tfo.writeln(Strings.concat(header, Strings.tab));
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			
			Integer topsnp = Integer.parseInt(elems[topsnpcol]);
			if (topsnp > 0) {
				tfo.writeln(Strings.concat(elems, Strings.tab));
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		
		tfo.close();
		tf.close();
	}
	
	public void run(String allelefile, String efile, String outefile) throws IOException {
		
		HashMap<String, String> refalleles = new HashMap<String, String>();
		HashMap<String, String> altalleles = new HashMap<String, String>();
		
		TextFile tf = new TextFile(allelefile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[0];
			String ref = elems[elems.length - 2];
			String alt = elems[elems.length - 1];
			
			refalleles.put(snp, ref);
			altalleles.put(snp, alt);
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile outf = new TextFile(outefile, TextFile.W);
		TextFile inf = new TextFile(efile, TextFile.R);
		outf.writeln(inf.readLine() + "\tRef\tAlt");
		elems = inf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[7];
			String ref = refalleles.get(snp);
			String alt = altalleles.get(snp);
			
			outf.writeln(Strings.concat(elems, Strings.tab) + "\t" + ref + "\t" + alt);
			elems = inf.readLineElems(TextFile.tab);
		}
		inf.close();
		outf.close();
		
	}
	
}
