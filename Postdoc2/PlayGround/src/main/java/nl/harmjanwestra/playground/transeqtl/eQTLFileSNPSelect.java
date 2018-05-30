package nl.harmjanwestra.playground.transeqtl;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class eQTLFileSNPSelect {
	
	public static void main(String[] args) {
		eQTLFileSNPSelect s = new eQTLFileSNPSelect();
		String[] infiles = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-01-31-cis-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-04-03-trans-eQTLsPrunedFDR0.05.txt.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\conditional_trans_eQTL_BIOSplusEGCUT_20180418-cond1.txt"
		};
		
		String[] outfilesplink = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-01-31-cis-snps-plink.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-04-03-trans-snps-plink.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\conditional_trans_eQTL_BIOSplusEGCUT_20180418-cond1-snps-plink.txt"
		};
		
		String[] outfilesproxyfinder = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-01-31-cis-snps-pf.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-04-03-trans-snps-pf.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\conditional_trans_eQTL_BIOSplusEGCUT_20180418-cond1-snps-pf.txt"
			
		};
		
		
		try {
			s.makesnplist(infiles, outfilesplink, outfilesproxyfinder);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		String[] prunedfiles = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-01-31-cis-snps-plink-pruned.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-04-03-trans-snps-plink-pruned.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-04-18-trans-cond-snps-plink-pruned.txt"
		};
		String[] outfilesproxyfinderpruned = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-01-31-cis-snps-pf-pruned.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-04-03-trans-snps-pf-pruned.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-17-freeze\\2018-04-18-trans-cond-snps-pf-pruned.txt"
		};
		
		try {
			s.converToProxyFinderFiles(prunedfiles, outfilesproxyfinder, outfilesproxyfinderpruned);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void converToProxyFinderFiles(String[] prunedfiles, String[] outfilesproxyfinder, String[] outfilesproxyfinderpruned) throws IOException {
		for (int i = 0; i < outfilesproxyfinder.length; i++) {
			TextFile in = new TextFile(outfilesproxyfinder[i], TextFile.R);
			HashMap<String, String> vartovar = new HashMap<>();
			String[] elems = in.readLineElems(TextFile.tab);
			while (elems != null) {
				String rs = elems[2];
				vartovar.put(rs, Strings.concat(elems, Strings.tab));
				elems = in.readLineElems(TextFile.tab);
			}
			in.close();
			
			TextFile in2 = new TextFile(prunedfiles[i], TextFile.R);
			TextFile out = new TextFile(outfilesproxyfinderpruned[i], TextFile.W);
			String rs = in2.readLine();
			while (rs != null) {
				out.writeln(vartovar.get(rs));
				rs = in2.readLine();
			}
			in2.close();
			out.close();
		}
	}
	
	public void makesnplist(String[] input, String[] outputplink, String[] outputpf) throws IOException {
		for (int i = 0; i < input.length; i++) {
			TextFile in = new TextFile(input[i], TextFile.R);
			in.readLine();
			String[] elems = in.readLineElems(TextFile.tab);
			TextFile out1 = new TextFile(outputplink[i], TextFile.W);
			TextFile out2 = new TextFile(outputpf[i], TextFile.W);
			HashSet<String> seen = new HashSet<String>();
			while (elems != null) {
				if (!seen.contains(elems[1])) {
					out1.writeln(elems[1]);
					
					out2.writeln(elems[2] + "\t" + elems[3] + "\t" + elems[1]);
					seen.add(elems[1]);
				}
				elems = in.readLineElems(TextFile.tab);
			}
			in.close();
			out1.close();
			out2.close();
		}
	}
}
