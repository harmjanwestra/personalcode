package nl.harmjanwestra.playground.trityper;

import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class CreateTestDataset {


	public static void main(String[] args) {

		String outdir = "d:\\tmp\\testdata\\";
		int numsamples = 333;
		CreateTestDataset c = new CreateTestDataset();
		try {
			c.run(outdir, numsamples);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String outdir, int numsamples) throws IOException {

		int[] skip = new int[]{5, 9, 11, 50, 75, 97, 153, 167, 168, 173, 190, 251};

		HashSet<Integer> skiphash = new HashSet<Integer>();
		for (int i : skip) {
			skiphash.add(i);
		}

		BinaryFile bf = new BinaryFile(outdir + "GenotypeMatrix.dat", BinaryFile.W);
		int nrpergt = numsamples / 3;

		TextFile sout = new TextFile(outdir + "Individuals.txt", TextFile.W);
		TextFile sout2 = new TextFile(outdir + "PhenotypeInformation.txt", TextFile.W);
		int indctr = 0;
		ArrayList<String> inds = new ArrayList<String>();

		for (int i = 0; i < nrpergt * 3; i++) {
			sout.writeln("Ind" + indctr);
			sout2.writeln("Ind" + indctr + "\tunknown\tinclude\tunknown");
			inds.add("Ind" + indctr);
			indctr++;
		}
		sout2.close();
		sout.close();


		// write genotypes
		byte[] a1 = new byte[nrpergt * 3];
		byte[] a2 = new byte[nrpergt * 3];

		byte al1 = 67;
		byte al2 = 69;

		for (int i = 0; i < nrpergt * 3; i++) {
			if (skiphash.contains(i)) {
				a1[i] = 0;
				a2[i] = 0;
			} else {
				if (i < nrpergt) {
					a1[i] = al1;
					a2[i] = al1;
				} else if (i < nrpergt * 2) {
					a1[i] = al1;
					a2[i] = al2;
				} else {
					a1[i] = al2;
					a2[i] = al2;
				}
			}
		}

		bf.write(a1);
		bf.write(a2);
		bf.close();

		TextFile snpout = new TextFile(outdir + "SNPs.txt", TextFile.W);
		TextFile snpout2 = new TextFile(outdir + "SNPMappings.txt", TextFile.W);
		snpout.writeln("rs1");
		snpout2.writeln("1\t1\trs1");
		snpout.close();
		snpout2.close();

		TextFile expout = new TextFile(outdir + "ExpressionData.txt", TextFile.W);

		String header = "-\t" + Strings.concat(inds, Strings.tab);
		expout.writeln(header);
		String gene = "ENSG01";
		for (int i = 0; i < inds.size(); i++) {
			gene += "\t" + (i + 1);
		}
		expout.writeln(gene);
		expout.close();

		TextFile annotout = new TextFile(outdir + "ExpressionAnnotation.txt", TextFile.W);
		header = "Platform\tArrayAddress\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tSeq";
		annotout.writeln(header);
		String outln = "test\tENSG01\tENSG01\t1\t1000\t1500\tENSG01\t-";
		annotout.writeln(outln);
		annotout.close();

	}
}
