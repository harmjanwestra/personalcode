package nl.harmjanwestra.playground.ase;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;

import java.io.IOException;
import java.util.HashMap;

public class CompareASEEffectDirectionWithASE {

	public static void main(String[] args) {

		String eqtl = "D:\\tmp\\niek\\eQTLsFDR-ProbeLevel.txt.gz";
		String ase = "D:\\tmp\\niek\\ase_and_eqtl.txt";
		String out = "D:\\tmp\\niek\\ase_and_eqtl_comp.txt";
		CompareASEEffectDirectionWithASE c = new CompareASEEffectDirectionWithASE();
		try {
			c.run(eqtl, ase, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public class EQTL {
		String snp;
		String gene;
		String alleles;
		String assessed;
		double z;
	}

	public void run(String eqtl, String ase, String outf) throws IOException {

		HashMap<String, EQTL> eqtls = new HashMap<>();

		TextFile tf = new TextFile(eqtl, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[2] + ":" + elems[3];
			String gene = elems[4];
			String alleles = elems[8];
			String assessed = elems[9];
			Double z = Double.parseDouble(elems[10]);

			EQTL e = new EQTL();
			e.snp = snp;
			e.alleles = alleles;
			e.assessed = assessed;
			e.z = z;
			e.gene = gene;
			eqtls.put(snp + "-" + gene, e);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile out = new TextFile(outf, TextFile.W);
		out.writeln("snp\tgene\tRatioASE\tAllelesASE\tAssessedASE\tZEQTL\tZEQTLFlip\tAllelesEQTL\tAssessedEQTL\tFlip");
		TextFile tf2 = new TextFile(ase, TextFile.R);
		tf2.readLine();
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[1] + ":" + elems[2];
			String gene = elems[6];
			String alleles = elems[3] + "/" + elems[4];
			String assessed = elems[4];
			Double z = Double.parseDouble(elems[11]);
			EQTL e = eqtls.get(snp + "-" + gene);
			if (e != null) {
				Boolean flip = BaseAnnot.flipalleles(alleles, assessed, e.alleles, e.assessed);
				Double zflip = e.z;
				if (flip == null) {
					zflip = Double.NaN;
				} else if (flip) {
					zflip = zflip * -1;
				}
				out.writeln(snp + "\t" + gene + "\t" + z + "\t" + alleles + "\t" + assessed + "\t" + e.z + "\t" + zflip + "\t" + e.alleles + "\t" + e.assessed + "\t" + flip);
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		out.close();
		tf2.close();

	}
}
