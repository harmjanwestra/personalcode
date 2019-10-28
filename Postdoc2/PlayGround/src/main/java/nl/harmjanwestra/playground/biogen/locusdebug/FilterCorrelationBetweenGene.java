package nl.harmjanwestra.playground.biogen.locusdebug;

import com.itextpdf.text.DocumentException;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class FilterCorrelationBetweenGene {

	public static void main(String[] args) {
		String eqtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Cis-Genes\\Primary\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
		String patchgenelist = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\2019-09-26-patch_genes.txt";
		String correlationfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\EUR-run2-dump-CorrelationBetweenGenes.txt";

		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\EUR-run2-dump-CorrelationBetweenGenes-significantgenes.pdf";
		FilterCorrelationBetweenGene c = new FilterCorrelationBetweenGene();
		try {
			c.run(eqtlfile, patchgenelist, correlationfile, output);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}

	}

	public class Gene {
		String name;
		HashSet<String> ensgids = new HashSet<String>();
		boolean hasautosomalcopies = false;
		boolean haspatchcopies = false;
	}

	public void run(String eqtlfile, String patchgenelistfile, String correlationfile, String output) throws IOException, DocumentException {


		HashSet<String> significantgenes = new HashSet<String>();

		TextFile tf = new TextFile(eqtlfile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			significantgenes.add(elems[4]);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		TextFile tf2 = new TextFile(patchgenelistfile, TextFile.R);
		String[] header = tf2.readLineElems(TextFile.tab);

		HashMap<String, Gene> ensgToGene = new HashMap<String, Gene>();
		ArrayList<Gene> genes = new ArrayList<>();
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {

			Gene g = new Gene();
			g.name = elems[0];

			for (int i = 1; i < elems.length; i++) {
				if (!elems[i].equals("-")) {
					String[] ensgids = elems[i].split(";");
					for (String s : ensgids) {
						g.ensgids.add(s);
						ensgToGene.put(s, g);
					}
					if (i < 23 && g.ensgids.size() > 1) {
						g.hasautosomalcopies = true;
					} else if (i >= 23) {
						g.haspatchcopies = true;
					}
				}
			}
			genes.add(g);

			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();

		ArrayList<ArrayList<Double>> xs = new ArrayList<>();
		ArrayList<ArrayList<Double>> ys = new ArrayList<>();
		for (int i = 0; i < 3; i++) {
			xs.add(new ArrayList<>());
			ys.add(new ArrayList<>());
		}


		TextFile tf3 = new TextFile(correlationfile, TextFile.R);
		tf3.readLine();
		elems = tf3.readLineElems(TextFile.tab);
		while (elems != null) {
			String gene = elems[0];
			if (significantgenes.contains(gene)) {
				Gene g = ensgToGene.get(gene);
				double v1 = Double.parseDouble(elems[1]);
				double v2 = Double.parseDouble(elems[2]);
				if (g != null) {
					if (g.haspatchcopies) {
						xs.get(2).add(v1);
						ys.get(2).add(v2);
					} else {
						xs.get(1).add(v1);
						ys.get(1).add(v2);
					}
				} else {
					xs.get(0).add(v1);
					ys.get(0).add(v2);
				}
			}
			elems = tf3.readLineElems(TextFile.tab);
		}
		tf3.close();

		double[][] xarr = new double[3][];
		double[][] yarr = new double[3][];
		for (int i = 0; i < xarr.length; i++) {

			xarr[i] = Primitives.toPrimitiveArr(xs.get(i));
			yarr[i] = Primitives.toPrimitiveArr(ys.get(i));
		}
		String[] labels = new String[]{"Autosomal", "AutosomalWithCopies", "PatchCopies"};

		Grid grid = new Grid(500, 500, 1, 1, 100, 100);
		ScatterplotPanel p = new ScatterplotPanel(1, 1);
		p.setData(xarr, yarr);
		p.setDatasetLabels(labels);
		p.setPlotElems(true, true);
		p.setDataRange(new Range(-1, -1, 1, 1));
		grid.addPanel(p);
		grid.draw(output);

	}

}
