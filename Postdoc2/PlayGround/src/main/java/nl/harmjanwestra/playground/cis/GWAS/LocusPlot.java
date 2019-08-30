package nl.harmjanwestra.playground.cis.GWAS;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.features.Gene;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.panels.AssociationPanel;
import umcg.genetica.graphics.panels.GenePanel;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.util.Primitives;
import umcg.genetica.util.RunTimer;
import umontreal.iro.lecuyer.util.Num;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

public class LocusPlot {


	public static void main(String[] args) {

		String genefile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\eqtlmerge\\MS\\ENSG00000013725.txt.gz";
		String colocfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\eqtlmerge\\MS\\ENSG00000013725.txt.gz-coloc.gz";
		String gtf = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\annotation\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz";
		String outputplot = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\eqtlmerge\\ms-ENSG00000013725-coloc.pdf";
		LocusPlot lp = new LocusPlot();
		try {
			lp.run(genefile, colocfile, gtf, outputplot);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
		RunTimer
	}


	public void run(String geneFile, String colocfile, String gtf, String outputplot) throws IOException, DocumentException {

		int maxpos = 0;
		int minpos = Integer.MAX_VALUE;
		Chromosome chr = null;
		HashMap<String, Double> colocresults = readColocFile(colocfile);

		ArrayList<Double> positions = new ArrayList<Double>();
		ArrayList<Double> z1s = new ArrayList<>();
		ArrayList<Double> z2s = new ArrayList<>();
		ArrayList<Double> pp4 = new ArrayList<>();

		TextFile tf = new TextFile(geneFile, TextFile.R);
		tf.readLine();

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[0];
			Double d = colocresults.get(snp);
			if (d != null) {
				String[] snpelems = snp.split(":");
				chr = Chromosome.parseChr(snpelems[0]);
				int pos = Integer.parseInt(snpelems[1]);
				if (pos > maxpos) {
					maxpos = pos;
				}
				if (pos < minpos) {
					minpos = pos;
				}

				positions.add((double) pos);
				double b1 = Double.parseDouble(elems[1]);
				double se1 = Double.parseDouble(elems[2]);
				double z1 = b1 / se1;
				z1s.add(z1);
				double b2 = Double.parseDouble(elems[5]);
				double se2 = Double.parseDouble(elems[6]);
				double z2 = b2 / se2;
				z2s.add(z2);
				pp4.add(d);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		GTFAnnotation annot = new GTFAnnotation(gtf);

		Collection<Gene> genes = annot.getGenes();
		Feature f = new Feature(chr, minpos, maxpos);
		ArrayList<Gene> genesInRegion = new ArrayList<>();
		for (Gene g : genes) {
			if (g.overlaps(f)) {
				genesInRegion.add(g);
			}
		}


		Grid g = new Grid(300, 200, 5, 1, 100, 100);

		GenePanel gp = new GenePanel(1, 1);
		gp.setData(f, genesInRegion);

		ScatterplotPanel pz1 = new ScatterplotPanel(1, 1);
		pz1.setData(Primitives.toPrimitiveArr(positions), Primitives.toPrimitiveArr(z1s));
		pz1.setTitle("MetaBrain cortex eQTL");
		pz1.setAlpha(0.2f);
		pz1.setLabels("Position (bp)", "eQTL Z-score");
		pz1.setPlotElems(true, false);

		ScatterplotPanel pz2 = new ScatterplotPanel(1, 1);
		pz2.setData(Primitives.toPrimitiveArr(positions), Primitives.toPrimitiveArr(z2s));
		pz2.setTitle("MS GWAS");
		pz2.setAlpha(0.2f);
		pz2.setLabels("Position (bp)", "GWAS Z-score");
		pz2.setPlotElems(true, false);

		ScatterplotPanel ppp4 = new ScatterplotPanel(1, 1);
		ppp4.setData(Primitives.toPrimitiveArr(positions), Primitives.toPrimitiveArr(pp4));
		ppp4.setTitle("COLOC colocalization");
		ppp4.setLabels("Position (bp)", "COLOC P(h4)");
		ppp4.setAlpha(0.8f);
		ppp4.setPlotElems(true, false);

		g.addPanel(gp);
		g.addPanel(pz1);
		g.addPanel(pz2);
		g.addPanel(ppp4);

		g.draw(outputplot);


	}

	private HashMap<String, Double> readColocFile(String colocfile) throws IOException {
		TextFile tf = new TextFile(colocfile, TextFile.R);
		tf.readLine();
		HashMap<String, Double> output = new HashMap<>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			try {
				String snp = elems[1];
				Double d = Double.parseDouble(elems[elems.length - 1]);
				output.put(snp, d);
			} catch (NumberFormatException e) {

			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();
		return output;
	}

}
