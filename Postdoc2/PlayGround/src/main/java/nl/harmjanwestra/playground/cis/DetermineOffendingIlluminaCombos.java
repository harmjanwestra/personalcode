package nl.harmjanwestra.playground.cis;

import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class DetermineOffendingIlluminaCombos {
	
	
	public static void main(String[] args) {
		
		if (args.length < 4) {
			System.out.println("Usage: illuminadatasets eqtlfile allowedcombos out");
			System.exit(-1);
		}
		
		String illuminaDatasetIds = args[0];
		String eqtlfile = args[1];
		String allowedcombos = args[2];
		String outf = args[3];
		
		DetermineOffendingIlluminaCombos c = new DetermineOffendingIlluminaCombos();
		
		try {
			c.run(illuminaDatasetIds, eqtlfile, allowedcombos, outf);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String illuminaDatasetIDs, String eQTLFile, String allowedCombos, String outf) throws IOException {
		
		THashMap<String, THashSet<String>> combos = new THashMap<>();
		TextFile tf = new TextFile(allowedCombos, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ct = 0;
		while (elems != null) {
			String snp = Strings.cache(elems[0]);
			String gene = Strings.cache(elems[1]);
			THashSet<String> hashset = combos.get(gene);
			if (hashset == null) {
				hashset = new THashSet<>();
			}
			hashset.add(snp);
			combos.put(gene,hashset);
			
			ct++;
			if (ct % 100000 == 0) {
				System.out.print(ct + " combos read.\r");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println(combos.size() + " combos allowed");
		
		TextFile tf2 = new TextFile(illuminaDatasetIDs, TextFile.R);
		HashSet<String> illuminads = new HashSet<String>();
		String ln = tf2.readLine();
		while (ln != null) {
			illuminads.add(ln);
			ln = tf2.readLine();
		}
		tf2.close();
		
		System.out.println(illuminads.size() + " illumina datasets");
		
		
		TextFile out = new TextFile(outf, TextFile.W);
		TextFile in = new TextFile(eQTLFile, TextFile.R);
		out.writeln(in.readLine());
		ln = in.readLine();
		int ctr = 0;
		int ctr2 = 0;
		int sign = 0;
		while (ln != null) {
			elems = Strings.tab.split(ln);
			
			
			
			String gene = elems[4];
			THashSet<String> snps = combos.get(gene);
			if (snps!=null && !snps.contains(elems[1])) {
				EQTL e = EQTL.fromString(elems, "-", Strings.semicolon);
				String[] dataset = e.getDatasets();
				
				boolean hasIllumina = false;
				for (String ds : dataset) {
					if (illuminads.contains(ds)) {
						hasIllumina = true;
						break;
					}
				}
				
				if (hasIllumina) {
					out.writeln(ln);
					ctr2++;
					if (e.getFDR() < 0.05) {
						sign++;
					}
				}
			}
			ln = in.readLine();
			ctr++;
			if (ctr % 100000 == 0) {
				System.out.println(ctr + " lines processed. " + ctr2 + " lines should possibly not be there, of which " + sign + " are significant.");
			}
		}
		
		System.out.println("done");
		System.out.println();
		System.out.println(ctr + " lines processed. " + ctr2 + " lines should possibly not be there, of which " + sign + " are significant.");
		in.close();
		out.close();
	}
}
