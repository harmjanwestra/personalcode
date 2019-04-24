package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class FilterGWASCatalog {

	public static void main(String[] args) {
		String gwasin = "D:\\Sync\\SyncThing\\Data\\Ref\\gwascatalog\\gwas_catalog_v1.0-associations_e95_r2019-03-22.tsv.gz";
		String gwasout = "D:\\Sync\\SyncThing\\Data\\Ref\\gwascatalog\\gwas_catalog_v1.0-associations_e95_r2019-03-22-snps.txt.gz";
		String gwasoutdisease = "D:\\Sync\\SyncThing\\Data\\Ref\\gwascatalog\\gwas_catalog_v1.0-associations_e95_r2019-03-22-snpsAndDiseases.txt.gz";

		try {
			GWASCatalog catalog = new GWASCatalog(gwasin);

			HashSet<GWASSNP> snps = catalog.getSnps();

			TextFile out1 = new TextFile(gwasout, TextFile.W);
			TextFile out2 = new TextFile(gwasoutdisease, TextFile.W);
			for (GWASSNP snp : snps) {
				out1.writeln(snp.getName());
				GWASTrait[] traits = snp.getAssociatedTraitsArray();
				String[] diseaseStr = new String[traits.length];
				for (int i = 0; i < diseaseStr.length; i++) {
					diseaseStr[i] = traits[i].getCleanName();
				}

				out2.writeln(snp.getName() + "\t" + Strings.concat(diseaseStr, Strings.semicolon));
			}
			out1.close();
			out2.close();

		} catch (IOException e) {
			e.printStackTrace();
		}


	}
}
