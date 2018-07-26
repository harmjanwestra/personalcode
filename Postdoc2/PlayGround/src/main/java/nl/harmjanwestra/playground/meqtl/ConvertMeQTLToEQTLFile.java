package nl.harmjanwestra.playground.meqtl;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class ConvertMeQTLToEQTLFile {
	
	
	public static void main(String[] args) {
		String eqtl = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\trans\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05.txt.gz";
		String eqtm = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-04-10-Replication\\2018-04-10-cis-eqtm-1mb.txt";
		String snpcgfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-04-10-Replication\\2018-07-25-doublecheck\\snpcgpairs.txt";
		
		String meqtl = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-04-10-Replication\\2017-04-10-trans-meqtl-transrepl.txt.gz";
		String replicationout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-04-10-Replication\\2018-07-25-doublecheck\\transeQTLToMeQTL.txt";
		
		ConvertMeQTLToEQTLFile m = new ConvertMeQTLToEQTLFile();
		
		try {
			m.eqtmtoeqtl(eqtm, eqtl, snpcgfile);
			
			m.eqtlreplicate(eqtl, eqtm, meqtl, replicationout);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private void eqtlreplicate(String eqtl, String eqtm, String meqtl, String replicationout) throws IOException {
		
		TextFile tf = new TextFile(eqtm, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		
		HashMap<String, EQTL> geneToCG = new HashMap<String, EQTL>();
		while (elems != null) {
			String gene = elems[4];
			if (!geneToCG.containsKey(gene)) {
				EQTL e = EQTL.fromString(elems, "-", Strings.semicolon);
				geneToCG.put(gene, e);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println(geneToCG.size() + " cgs loaded");
		
		HashMap<String, EQTL> meqtls = new HashMap<>();
		
		TextFile tf2 = new TextFile(meqtl, TextFile.R);
		tf2.readLine();
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			EQTL e = EQTL.fromString(elems, "-", Strings.semicolon);
			String query = e.getRsName() + "_" + e.getProbe();
			if (!meqtls.containsKey(query)) {
				meqtls.put(query, e);
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		
		System.out.println(meqtls.size() + " meqtls loaded ");
		TextFile tf3 = new TextFile(eqtl, TextFile.R);
		TextFile out = new TextFile(replicationout, TextFile.W);
		String header = "Gene\tSNP\tAlleles\tAssessed\tZ\tP\tFDR\tEQTM-CG\tEQTM-Alleles\tEQTM-Assessed\tEQTM-Z\tEQTM-P\tEQTM-FDR\tMEQTL-Alleles\tMEQTL-Assessed\tMEQTL-Z\tMEQTL-P\tMEQTL-FDR";
		out.writeln(header);
		tf3.readLine();
		elems = tf3.readLineElems(TextFile.tab);
		int read = 0;
		int merged = 0;
		while (elems != null) {
			EQTL e = EQTL.fromString(elems, "-", Strings.semicolon);
			
			EQTL em = geneToCG.get(e.getProbe());
			EQTL me = null;
			if (em != null) {
				
				String query = e.getRsName() + "_" + em.getRsName();
				me = meqtls.get(query);
			}
			
			
			String outln = e.getProbe()
					+ "\t" + e.getRsName()
					+ "\t" + e.getAlleles()
					+ "\t" + e.getAlleleAssessed()
					+ "\t" + e.getZscore()
					+ "\t" + e.getPvalue()
					+ "\t" + e.getFDR();
			
			if (em == null) {
				outln += "\t-\t-\t-\t-\t-\t-";
				outln += "\t-\t-\t-\t-\t-\t-";
			} else {
				outln += "\t" + em.getRsName()
						+ "\t" + em.getAlleles()
						+ "\t" + em.getAlleleAssessed()
						+ "\t" + em.getZscore()
						+ "\t" + em.getPvalue()
						+ "\t" + em.getFDR();
				
				if (me == null) {
					outln += "\t-\t-\t-\t-\t-\t-";
				} else {
					outln += "\t" + me.getRsName()
							+ "\t" + me.getAlleles()
							+ "\t" + me.getAlleleAssessed()
							+ "\t" + me.getZscore()
							+ "\t" + me.getPvalue()
							+ "\t" + me.getFDR();
					merged++;
				}
			}
			out.writeln(outln);
			read++;
			elems = tf3.readLineElems(TextFile.tab);
		}
		out.close();
		tf3.close();
		
		System.out.println(read + " eqtls merged to " + merged + " meqtl");
		
	}
	
	public void eqtmtoeqtl(String eqtmfile, String eqtlfile, String output) throws IOException {
		
		TextFile tf = new TextFile(eqtmfile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		
		HashMap<String, String> geneToCG = new HashMap<String, String>();
		while (elems != null) {
			String cg = elems[1];
			String gene = elems[4];
			if (!geneToCG.containsKey(gene)) {
				geneToCG.put(gene, cg);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile tf2 = new TextFile(eqtlfile, TextFile.R);
		tf2.readLine();
		elems = tf2.readLineElems(TextFile.tab);
		TextFile outf = new TextFile(output, TextFile.W);
		while (elems != null) {
			String snp = elems[1];
			String gene = elems[4];
			String cg = geneToCG.get(gene);
			if (cg != null) {
				outf.writeln(snp + "\t" + cg + "\t" + gene);
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		outf.close();
		tf2.close();
		
		
	}
	
	
}
