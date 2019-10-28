package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class DiseaseEnrichmentComparison {


	public static void main(String[] args) {
		DiseaseEnrichmentComparison c = new DiseaseEnrichmentComparison();
		String gwascatalog = "D:\\Sync\\SyncThing\\Data\\Ref\\gwascatalog\\gwas_catalog_v1.0-associations_e96_r2019-10-14.tsv.gz";
		String snpmap = "D:\\Sync\\SyncThing\\Data\\Ref\\dbsnp\\SNPMappings-dbsnp151-b38.txt.gz";
		String gwascatalogpositions = "D:\\Sync\\SyncThing\\Data\\Ref\\gwascatalog\\gwas_catalog_v1.0-associations_e96_r2019-10-14-b38positions.txt.gz";


		try {
			/* procedure:
		       0. read gwas catalog
		       1. map snps to positions (b38)
		    */
			c.mapGWASCatalogToSNPPositions(gwascatalog, gwascatalogpositions, snpmap, 5e-8);


			String eqtlgengenelist = "";
			String metabraingenelist = "";
//			c.overlapGeneLists(eqtlgengenelist, metabraingenelist);



		} catch (IOException e) {
			e.printStackTrace();
		}

		/* procedure:
		0. read gwas catalog
		1. map snps to positions (b38)

		1. determine overlap in tested genes
		2. determine overlap in tested snps // is this required?

		for each disease:
		    1. determine snps associated with disease
		    2. determine set of 'other' snps
		    3. limit snps to those tested in both datasets // may not be required
		    for each eqtl dataset:
		        for each gene:
		            1. determine whether gene is tested in both datasets
		            2. determine cis-region
		            2. determine disease snps in cis-region
		            3. determine ld between eqtl
		      2.

		 */

	}

	public void mapGWASCatalogToSNPPositions(String catalog, String output, String snpmap, double threshold) throws IOException {
		GWASCatalog c = new GWASCatalog(catalog);

		HashSet<String> rsids = new HashSet<String>();

		GWASTrait[] traits = c.getTraits();

		for (GWASTrait t : traits) {
			GWASSNP[] snps = t.getSNPs(threshold);
			int inc = 0;
			int excl = 0;
			for (GWASSNP snp : snps) {
				if (snp.getName().contains("rs")) {
					rsids.add(snp.getName());
					inc++;
				} else {
					System.out.println("Chose to exclude: " + snp.getName());
					excl++;
				}
			}
			System.out.println(t.getName() + "\t" + inc + " included, " + excl + " excluded variants.");
		}

		System.out.println(rsids.size() + " rs ids found..");
		TextFile tf = new TextFile(snpmap, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, String> rsIdToPos = new HashMap<>();
		int ctr = 0;
		while (elems != null) {
			String rs = elems[2];
			if (rsids.contains(rs)) {
				rsIdToPos.put(rs, elems[0] + ":" + elems[1]);
			}
			ctr++;
			if (ctr % 1000000 == 0) {
				System.out.println(ctr + " snp rsids parsed; " + rsIdToPos.size() + " / " + rsids.size() + " mapped.");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println("Could recover " + rsIdToPos.size() + " / " + rsids.size() + " rs ids");

		TextFile out = new TextFile(output, TextFile.W);

		out.writeln("Trait\tnrSNPs\tSNPs");
		for (GWASTrait t : traits) {
			GWASSNP[] snps = t.getSNPs(threshold);
			String ln = t.getName();
			ArrayList<String> positions = new ArrayList<String>();
			for (GWASSNP snp : snps) {
				if (snp.getName().contains("rs")) {
					String pos = rsIdToPos.get(snp.getName());
					if (pos != null) {
						pos = pos + ":" + snp.getName();
						positions.add(pos);
					}
				}
			}
			if (!positions.isEmpty()) {
				out.writeln(ln + "\t" + positions.size() + "\t" + Strings.concat(positions, Strings.semicolon));
			}
		}
		out.close();


	}


	public void run() {


	}
}
