package nl.harmjanwestra.playground.eqtlgen;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;

public class CompareBloodWithPBMC {

	public static void main(String[] args) {
		String in = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-09-03-PBMCvsWholeBlood\\2018-04-10-transEQTL-Replication-PBMCAndWholeBlood.txt.gz";
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-09-03-PBMCvsWholeBlood\\2018-04-10-transEQTL-Replication-PBMCAndWholeBlood-groups.txt";
		CompareBloodWithPBMC p = new CompareBloodWithPBMC();
		try {
			p.run(in, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String filein, String fileout) throws IOException {
		TextFile f = new TextFile(filein, TextFile.R);

		f.readLine(); // header 1
		f.readLine(); // header 2

		int genecol = 0;
		int snpcol = 3;
		int pbmczcol = 12;
		int pbmcfdrcol = 16;
		int wbzcol = 17;
		int wbfdrcol = 21;


		String[] elems = f.readLineElems(TextFile.tab);
		TextFile outf = new TextFile(fileout, TextFile.W);
		outf.writeln("SNP\tGene\tCategory");
		while (elems != null) {

			String gene = elems[genecol];
			String snp = elems[snpcol];
			String pbmczs = elems[pbmczcol];
			String pbmcfs = elems[pbmcfdrcol];
			String wbzs = elems[wbzcol];
			String wbfs = elems[wbfdrcol];

			double pbmcz = 0;
			double pbmcf = 2;
			double wbz = 0;
			double wbf = 2;
			try {
				pbmcz = Double.parseDouble(pbmczs);
				pbmcf = Double.parseDouble(pbmcfs);
				wbz = Double.parseDouble(wbzs);
				wbf = Double.parseDouble(wbfs);
			} catch (NumberFormatException e) {

			}


			// PBMC: no ery's etc
			// whole blood: all blood cells
			// effects present in wb but not pbmc are likely driven by cell type
			// effects present in pbmc but not wb are less likely driven by cell type
			// shared effects are less likely to be cell driven

			if (wbf > 1 || pbmcf > 1) {
				outf.writeln(snp + "\t" + gene + "\tNotPresent");
			} else {
				if (pbmcf < 0.05 && wbf < 0.05) {
					// significant in both
					outf.writeln(snp + "\t" + gene + "\tPBMCandWholeBlood");
				} else if (pbmcf < 0.05 && wbf >= 0.05) {
					// significant in pbmc
					outf.writeln(snp + "\t" + gene + "\tPBMC");
				} else if (pbmcf >= 0.05 && wbf < 0.05) {
					// significant in wb
					outf.writeln(snp + "\t" + gene + "\tWholeBlood");
				} else {
					outf.writeln(snp + "\t" + gene + "\tNeither");
				}
			}

			elems = f.readLineElems(TextFile.tab);
		}
		f.close();
		outf.close();

	}
}
