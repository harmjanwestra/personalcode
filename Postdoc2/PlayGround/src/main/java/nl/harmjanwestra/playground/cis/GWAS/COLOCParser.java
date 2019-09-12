package nl.harmjanwestra.playground.cis.GWAS;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;

public class COLOCParser {

	public static void main(String[] aregs) {

		String dir = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\eqtlmerge\\ALZ\\";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\eqtlmerge\\ALZ-genes.txt";
		String output2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\eqtlmerge\\ALZ-genesWithColocSNPs.txt";
		COLOCParser p = new COLOCParser();
		try {
			p.parseSummariesInDir(dir, output);
			p.identifyLociWithColocalizingSNPs(dir, output2, 0.2);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void identifyLociWithColocalizingSNPs(String indir, String output, double threshold) throws IOException {
		String[] filelist = Gpio.getListOfFiles(indir);
		TextFile outf = new TextFile(output, TextFile.W);
		outf.writeln("Gene\tSNP\tPPH4");
		for (String filename : filelist) {
			if (filename.endsWith("coloc.gz")) {
				File f = new File(filename);
				String name = f.getName();
				String gene = name.replace(".txt.gz-coloc.gz", "");
				System.out.println("Parsing: " + filename);
				TextFile tf = new TextFile(indir + filename, TextFile.R);
				tf.readLine();
				String[] elems = tf.readLineElems(TextFile.tab);
				while (elems != null) {
					try {
						Double pp = Double.parseDouble(elems[elems.length - 1]);
						String snp = elems[1];
						if (pp > threshold) {
							outf.writeln(gene + "\t" + snp + "\t" + pp);
						}
					} catch (NumberFormatException e) {

					}
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
			}
		}
		outf.close();
	}

	public void parseSummariesInDir(String indir, String output) throws IOException {
		String[] filelist = Gpio.getListOfFiles(indir);
		TextFile outf = new TextFile(output, TextFile.W);
		outf.writeln("Gene\tnsnps\tPPH0\tPPH1\tPPH2\tPPH3\tPPH4");
		for (String filename : filelist) {
			if (filename.endsWith("colocsummary.gz")) {
				System.out.println("Parsing: " + filename);
				TextFile tf = new TextFile(indir + filename, TextFile.R);
				String ln = tf.readLine();
				String nsnps = "";
				String h0 = "";
				String h1 = "";
				String h2 = "";
				String h3 = "";
				String h4 = "";

				while (ln != null) {
					String[] elems = ln.split("\t");
					if (ln.startsWith("nsnps")) {
						nsnps = elems[1];
					} else if (ln.startsWith("PP.H0.abf")) {
						h0 = elems[1];
					} else if (ln.startsWith("PP.H1.abf")) {
						h1 = elems[1];
					} else if (ln.startsWith("PP.H2.abf")) {
						h2 = elems[1];
					} else if (ln.startsWith("PP.H3.abf")) {
						h3 = elems[1];
					} else if (ln.startsWith("PP.H4.abf")) {
						h4 = elems[1];
					}
					ln = tf.readLine();
				}
				tf.close();
				File f = new File(filename);
				String name = f.getName();
				String gene = name.replace(".txt.gz-colocsummary.gz", "");
				outf.writeln(gene + "\t" + nsnps + "\t" + h0 + "\t" + h1 + "\t" + h2 + "\t" + h3 + "\t" + h4);
			}
		}
		outf.close();
	}
}
