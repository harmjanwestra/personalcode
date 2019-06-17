package nl.harmjanwestra.playground.conv;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class TriTyperSNPFilter {


	public static void main(String[] args) {

		TriTyperSNPFilter s = new TriTyperSNPFilter();


		if (args.length < 3) {
			System.out.println("Usage: input snplist output [snploc] [snpmaploc]");
		} else {
			String snploc = null;
			String snpmaploc = null;
			if (args.length > 3) {
				snploc = args[3];
			}
			if (args.length > 4) {
				snpmaploc = args[4];
			}
			try {
				s.filter(args[0], args[1], MODE.INCLUDE, snploc, snpmaploc, args[2]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}

	enum MODE {
		INCLUDE,
		EXCLUDE
	}

	public void filter(String input, String snplist, MODE filtermode, String snploc, String snpmaploc, String output) throws IOException {

		if (snploc == null) {
			snploc = input + "SNPs.txt.gz";
		}
		if (snpmaploc == null) {
			snpmaploc = input + "SNPMappingss.txt.gz";
		}

		HashSet<String> querylist = new HashSet<String>();

		TextFile li = new TextFile(snplist, TextFile.R);
		querylist.addAll(li.readAsArrayList());
		li.close();

		System.out.println("Looking for " + querylist.size() + " snps ");

		Gpio.createDir(output);
		Gpio.copyFile(input + "Individuals.txt", output + "Individuals.txt");
		TextFile tfind = new TextFile(output + "Individuals.txt", TextFile.R);
		ArrayList<String> individuals = tfind.readAsArrayList();
		tfind.close();
		int nrInd = individuals.size();

		Gpio.copyFile(input + "PhenotypeInformation.txt", output + "PhenotypeInformation.txt");


		TextFile snpin = new TextFile(snploc, TextFile.R);
		TextFile snpout = new TextFile(output + "SNPs.txt.gz", TextFile.W);

		BinaryFile bfgti = new BinaryFile(input + "GenotypeMatrix.dat", BinaryFile.R);
		BinaryFile bfgto = new BinaryFile(output + "GenotypeMatrix.dat", BinaryFile.W);

		BinaryFile bfdosi = null;
		BinaryFile bfdoso = null;
		if (Gpio.exists(input + "ImputedDosageMatrix.dat")) {
			bfdosi = new BinaryFile(input + "ImputedDosageMatrix.dat", BinaryFile.R);
			bfdoso = new BinaryFile(output + "ImputedDosageMatrix.dat", BinaryFile.W);
		}


		String ln = snpin.readLine();
		byte[] gt = new byte[nrInd * 2];
		byte[] dos = new byte[nrInd];

		int ctr = 0;
		int written = 0;
		while (ln != null) {

			bfgti.read(gt);
			if (bfdosi != null) {
				bfdosi.read(dos);
			}
			if (querylist.contains(ln)) {
				if (filtermode.equals(MODE.INCLUDE)) {
					// write
					bfgto.write(gt);
					if (bfdosi != null) {
						bfdoso.write(dos);
					}
					snpout.writeln(ln);
					written++;
				}
			} else if (filtermode.equals(MODE.EXCLUDE)) {
				bfgto.write(gt);
				if (bfdosi != null) {
					bfdoso.write(dos);
				}
				snpout.writeln(ln);
				written++;
			}
			ln = snpin.readLine();

			ctr++;
			if (ctr % 100000 == 0) {
				System.out.print("\r" + written + " written \t" + ctr);
			}

		}
		System.out.println();
		System.out.println(written + "\t" + ctr);
		snpin.close();
		snpout.close();
		bfgti.close();
		bfgto.close();
		if (bfdosi != null) {
			bfdosi.close();
			bfdoso.close();
		}


		// filter snpmappings
		TextFile tfmapi = new TextFile(snpmaploc, TextFile.R);
		TextFile tfmapo = new TextFile(output + "SNPMappings.txt.gz", TextFile.W);
		String[] elems = tfmapi.readLineElems(TextFile.tab);
		while (elems != null) {

			String snps = elems[2];
			if (querylist.contains(snps)) {
				if (filtermode.equals(MODE.INCLUDE)) {
					tfmapo.writeln(Strings.concat(elems, Strings.tab));
				}
			} else if (filtermode.equals(MODE.EXCLUDE)) {
				tfmapo.writeln(Strings.concat(elems, Strings.tab));
			}
			elems = tfmapi.readLineElems(TextFile.tab);
		}

		tfmapo.close();
		tfmapi.close();
	}
}
