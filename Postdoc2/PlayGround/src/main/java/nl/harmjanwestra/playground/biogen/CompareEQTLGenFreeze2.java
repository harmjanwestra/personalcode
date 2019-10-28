package nl.harmjanwestra.playground.biogen;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.playground.legacy.vcf.DetermineLD;
import nl.harmjanwestra.playground.legacy.vcf.VCFTabix;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.SNPFeature;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class CompareEQTLGenFreeze2 {


	public static void main(String[] args) {

		String eqtlgen = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-03-EQTLGenRepl\\2018-01-31-eqtlgen-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38-topfxpergene.txt.gz";
		String metabrain = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-03-EQTLGenRepl\\eQTLsFDR-ProbeLevel.txt.gz";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-03-EQTLGenRepl\\output\\";

		String kg38 = "";
		String samplefilter = "";
		String eqtlgentop = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-03-EQTLGenRepl\\2018-01-31-eqtlgen-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38-topfxpergene.txt.gz";
		String metabraintop = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-03-EQTLGenRepl\\metabraintopfx\\EUR-run2-topfxPerGene.txt";

		CompareEQTLGenFreeze2 c = new CompareEQTLGenFreeze2();
		try {
			c.run(eqtlgen, metabrain, output);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}

	}

	public void ldmeh(String eqtlgen, String metabrain, String tabixmeh, String samplefilter, String output) throws IOException {

		HashMap<String, String> eqtlgensnps = new HashMap<String, String>();
		TextFile tf = new TextFile(eqtlgen, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			String snp = elems[2] + ":" + elems[3];
			String gene = elems[4];
			eqtlgensnps.put(gene, snp);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		TextFile tf2 = new TextFile(metabrain, TextFile.R);
		tf2.readLine();
		elems = tf2.readLineElems(TextFile.tab);

		TextFile outf = new TextFile(output, TextFile.W);
		outf.writeln("Gene\teqtlgensnp\tmetabrainsnp\trsq\tdprime");

		while (elems != null) {
			String snpchr = elems[2];
			String snppos = elems[3];
			String gene = elems[4];
			String countersnp = eqtlgensnps.get(gene);
			if (countersnp != null) {

				if (countersnp.equals(snpchr + ":" + snppos)) {

					// equal snp!
					outf.writeln(elems[4] + "\t" + elems[2] + ":" + elems[3] + "\t" + elems[2] + ":" + elems[3] + "\t" + 1d + "\t" + 1d);
				} else {
					String chrfile = tabixmeh.replaceAll("CHR", "" + snpchr);
					VCFTabix tabix = new VCFTabix(chrfile);
					boolean[] indfilter = tabix.getSampleFilter(samplefilter);

					int snppos1 = Integer.parseInt(snppos);
					SNPFeature snp = new SNPFeature(Chromosome.parseChr(snpchr), snppos1, snppos1);
					VCFVariant var1 = tabix.getVariant(snp, indfilter);

					String[] snp2elems = countersnp.split(":");

					int snppos2 = Integer.parseInt(snp2elems[1]);
					SNPFeature snp2 = new SNPFeature(Chromosome.parseChr(snp2elems[0]), snppos1, snppos1);
					VCFVariant var2 = tabix.getVariant(snp2, indfilter);

					DetermineLD ld = new DetermineLD();
					Pair<Double, Double> ldoutput = ld.getLD(var1, var2);

					double rsq = ldoutput.getRight();
					double dpr = ldoutput.getLeft();


					outf.writeln(elems[4] + "\t" + snp2elems[0] + ":" + snp2elems[1] + "\t" + elems[2] + ":" + elems[3] + "\t" + rsq + "\t" + dpr);

				}

			}


			elems = tf2.readLineElems(TextFile.tab);
		}


		tf2.close();


	}

	private class EQTL {
		double z;
		String allele;
		String assessed;
	}

	public void run(String eqtlgen, String metabrain, String output) throws IOException, DocumentException {
		HashMap<String, EQTL> eqtls = new HashMap<String, EQTL>();

		TextFile tf = new TextFile(eqtlgen, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			String id = elems[2] + ":" + elems[3] + ":" + elems[4];
			Double z = Double.parseDouble(elems[10]);
			String al = elems[8];
			String as = elems[9];
			EQTL eq = new EQTL();
			eq.z = z;
			eq.allele = al;
			eq.assessed = as;
			eqtls.put(id, eq);


			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		tf = new TextFile(metabrain, TextFile.R);
		tf.readLine();
		elems = tf.readLineElems(TextFile.tab);

//		ArrayList<Double> x = new ArrayList<Double>();
//		ArrayList<Double> y = new ArrayList<Double>();

		int overlap = 0;
		int overlapsig = 0;

		int eqdir = 0;
		int eqdirsig = 0;
		ScatterplotPanel p = new ScatterplotPanel(1, 1);
		ScatterplotPanel p2 = new ScatterplotPanel(1, 1);
		TextFile outf = new TextFile(output + "ZScoreCompEQTLgenVsMetaBrain.txt", TextFile.W);
		outf.writeln("EQTL\teQTLgenZ\tMetabrainZ\tConcordant\tSignificantMetabrain");
		while (elems != null) {


			String id = elems[2] + ":" + elems[3] + ":" + elems[4];
			Double z = Double.parseDouble(elems[10]);
			String al = elems[8];
			String as = elems[9];
			Double fdr = Double.parseDouble(elems[elems.length - 1]);

			EQTL eq = eqtls.get(id);

			if (eq != null) {

				Boolean flip = BaseAnnot.flipalleles(al, as, eq.allele, eq.assessed);
				if (flip != null) {
					overlap++;
					if (fdr < 0.05) {
						overlapsig++;
					}
					if (flip) {
						z *= -1;
					}

					boolean concordant = false;
					if ((z >= 0 && eq.z >= 0) || (z < 0 && eq.z < 0)) {
						eqdir++;
						concordant = true;

						if (fdr < 0.05) {
							eqdirsig++;
						}
					}
					boolean significant = false;
					if (fdr < 0.05) {
						significant = true;
					}
					p.addData(eq.z, z);
					outf.writeln(id + "\t" + z + "\t" + eq.z + "\t" + concordant + "\t" + significant);
					if (fdr < 0.05) {
						p2.addData(eq.z, z);
					}
				} else {
					System.out.println("Incompatible alleles for: " + id + "\t" + al + "\t" + as + "\t" + eq.allele + "\t" + eq.assessed);

				}
			} else {
				System.out.println("Could not find: " + id);
			}


			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		outf.close();

		System.out.println("Overlap: " + overlap);
		System.out.println("Overlap sig: " + overlapsig);
		System.out.println("EqDir: " + eqdir + " --> " + ((double) eqdir / overlap));
		System.out.println("EqDir sig: " + eqdirsig + " --> " + ((double) eqdirsig / overlapsig));

		Grid g = new Grid(500, 500, 1, 2, 100, 100);
		Range r = new Range(-40, -40, 40, 40);

		p.setDataRange(r);
		p2.setDataRange(r);
		p.setPlotElems(true, false);
		p2.setPlotElems(true, false);


		g.addPanel(p);
		g.addPanel(p2);

		g.draw(output + "ZScoreCompPlot.pdf");

	}
}
