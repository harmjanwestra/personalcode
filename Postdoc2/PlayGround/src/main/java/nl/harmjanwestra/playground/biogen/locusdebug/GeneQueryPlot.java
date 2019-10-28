package nl.harmjanwestra.playground.biogen.locusdebug;

import com.itextpdf.text.DocumentException;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.ViolinBoxPlot;
import umcg.genetica.graphics.panels.BoxPlotPanel;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.graphics.panels.SpacerPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.*;

public class GeneQueryPlot {

	public static void main(String[] args) {
		GeneQueryPlot p = new GeneQueryPlot();
		String inputmetabrain = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\MAPT-metabrain-Cis-Cortex-EUR.txt";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\MAPT-metabrain-Cis-Cortex-EUR-z.pdf";

		String[] gtexinput = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Amygdala-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Anterior_cingulate_cortex_BA24-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Caudate_basal_ganglia-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Cerebellar_Hemisphere-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Cerebellum-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Cortex-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Frontal_Cortex_BA9-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Hippocampus-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Hypothalamus-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Nucleus_accumbens_basal_ganglia-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Putamen_basal_ganglia-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Spinal_cord_cervical_c-1-ENSG00000186868.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\Brain_Substantia_nigra-ENSG00000186868.txt",
		};

		String[] gtexnames = new String[]{
				"Brain_Amygdala",
				"Brain_Anterior_cingulate_cortex_BA24",
				"Brain_Caudate_basal_ganglia",
				"Brain_Cerebellar_Hemisphere",
				"Brain_Cerebellum",
				"Brain_Cortex",
				"Brain_Frontal_Cortex_BA9",
				"Brain_Hippocampus",
				"Brain_Hypothalamus",
				"Brain_Nucleus_accumbens_basal_ganglia",
				"Brain_Putamen_basal_ganglia",
				"Brain_Spinal_cord_cervical_c-1",
				"Brain_Substantia_nigra"
		};

		String outputgtex = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\MAPT-metabrain-Cis-Cortex-EUR-GTEx.pdf";

		String expressiondata = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\exp\\20PCNoQQ-MAPT.txt";
		expressiondata = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\exp\\MAPT-NoDummyVars.txt";
		String sampleToDataset = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\GTE-EUR\\samplePerDataset.meh";
		String expoutput = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\exp\\20PCNoQQ-MAPT.pdf";
		 expoutput = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\MAPT\\exp\\MAPT-NoDummyVars.pdf";
		try {
			boolean zscores = false;
//			p.plotIndividualDatasetsEQTLGenFormat(inputmetabrain, output, zscores);
//			p.plotIndividualDatasetsGTExFormat(inputmetabrain, gtexinput, gtexnames, outputgtex);
			p.plotGeneExpression(expressiondata, sampleToDataset, "ENSG00000186868.15", expoutput);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void plotGeneExpression(String input, String samplesPerDataset, String gene, String output) throws Exception {

		HashSet<String> querygenes = new HashSet<String>();
		querygenes.add(gene);

		HashMap<String, String> sampleToDataset = new HashMap<>();
		HashSet<String> datasets = new HashSet<>();

		TextFile tf = new TextFile(samplesPerDataset, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, HashSet<String>> datasetToSample = new HashMap<>();
		while (elems != null) {
			sampleToDataset.put(elems[0], elems[1]);
			HashSet<String> samplesInDs = datasetToSample.get(elems[1]);
			if (samplesInDs == null) {
				samplesInDs = new HashSet<>();
			}
			samplesInDs.add(elems[0]);
			datasetToSample.put(elems[1], samplesInDs);
			datasets.add(elems[1]);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		HashSet<String> samples = new HashSet<String>();
		samples.addAll(sampleToDataset.keySet());

		ArrayList<String> dsList = new ArrayList<>();
		dsList.addAll(datasets);
		Collections.sort(dsList);

		DoubleMatrixDataset<String, String> data = DoubleMatrixDataset.loadSubsetOfTextDoubleData(input, '\t', querygenes, samples);

//		ViolinBoxPlot vbp = new ViolinBoxPlot();
		ViolinBoxPlot vbp = new ViolinBoxPlot();

		double[][][] dsVals = new double[dsList.size()][1][0];
		double[][] dsVals2 = new double[dsList.size()][0];
		String[][] xlabels = new String[dsList.size()][1];

		for (int d = 0; d < dsList.size(); d++) {
			xlabels[d][0] = "ds" + d;
			String ds = dsList.get(d);
			ArrayList<Double> vals = new ArrayList<>();
			HashSet<String> samplesInDs = datasetToSample.get(ds);
			for (String sample : samplesInDs) {
				Integer id = data.getHashCols().get(sample);
				if (id != null && id > -1) {
					vals.add(data.getElementQuick(0, id));
				}
			}
			dsVals[d][0] = Primitives.toPrimitiveArr(vals);
			dsVals2[d] = dsVals[d][0];
		}

		vbp.draw(dsVals, dsList.toArray(new String[0]), xlabels, gene, ViolinBoxPlot.Output.PDF, output);

		BoxPlotPanel p = new BoxPlotPanel(1, 1);
		p.setBinLabels(dsList.toArray(new String[0]));
		p.setData(dsVals2);
		p.setUseMeanAndSd(false);

		Grid g = new Grid(800, 300, 1, 1, 100, 100);
		g.addPanel(p);
		g.draw(output + "-grid.pdf");

	}

	class Dataset {
		ArrayList<Double> x = new ArrayList<>();
		ArrayList<Double> y = new ArrayList<>();
		ArrayList<Integer> n = new ArrayList<>();
		ArrayList<String> snps = new ArrayList<>();
		String name;

		double getMaxY() {
			double max = 0;
			for (double d : y) {
				if (d > max) {
					max = d;
				}
			}
			return max;
		}
	}

	public void plotIndividualDatasetsEQTLGenFormat(String input, String output, boolean useZ) throws IOException, DocumentException {


		ArrayList<Dataset> datasets = new ArrayList<Dataset>();
		TextFile tf = new TextFile(input, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		int maxpos = -1;
		int minpos = Integer.MAX_VALUE;
		Dataset meta = new Dataset();
		meta.name = "MetaAnalysis";
		ArrayList<String> snps = new ArrayList<>();
		while (elems != null) {

			if (elems.length > 1) {
				double metap = Double.parseDouble(elems[0]);
				double metaz = Double.parseDouble(elems[10]);
				String snp = elems[1];
				snps.add(snp);
				int snpchr = Integer.parseInt(elems[2]);
				int snppos = Integer.parseInt(elems[3]);
				if (snppos > maxpos) {
					maxpos = snppos;
				}
				if (snppos < minpos) {
					minpos = snppos;
				}
//				double metaz = Double.parseDouble(elems[10]);

				String[] ds = elems[11].split(";");

				String[] zs = elems[12].split(";");
				String[] ns = elems[13].split(";");
				int sum = 0;
				for (int i = 0; i < ds.length; i++) {
					if (ctr == 0) {
						datasets.add(new Dataset());
					}
					if (!ds[i].equals("-")) {
						String dsname = ds[i];
						Dataset dsobj = datasets.get(i);
						dsobj.name = dsname;
						dsobj.x.add((double) snppos);
						double dsz = Double.parseDouble(zs[i]);
						double dsp = -Math.log10(ZScores.zToP(dsz));
						if (useZ) {
							dsobj.y.add(dsz);
						} else {
							dsobj.y.add(dsp);
						}
						int n = Integer.parseInt(ns[i]);
						dsobj.n.add(n);
						sum += n;
					} else {
						Dataset dsobj = datasets.get(i);

						dsobj.x.add((double) snppos);

						dsobj.y.add(0d);

					}
				}

				meta.x.add((double) snppos);
				if (useZ) {

					meta.y.add(metaz);
				} else {
					double dsp = -Math.log10(metap);
					meta.y.add(dsp);
				}

				meta.n.add(sum);
			}

			ctr++;
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		int nrcol = 2;
		int nrrows = (datasets.size()) / nrcol;

		Grid grid = new Grid(500, 300, nrrows, 3, 100, 100);

		ScatterplotPanel p = new ScatterplotPanel(1, 1);
		p.setData(Primitives.toPrimitiveArr(meta.x), Primitives.toPrimitiveArr(meta.y));
		Range range = new Range(minpos, -40, maxpos, 40);
		range.round();
		range.setMinY(0d);
		range.setMinX(44922000);
		p.setDataRange(range);
		p.setTitle(meta.name);
		p.setPlotElems(true, false);
		grid.addPanel(p, 0, 0);
		SpacerPanel spacer = new SpacerPanel(1, 1);
		for (int r = 1; r < nrrows; r++) {
			grid.addPanel(spacer, r, 0);
		}

		for (int i = 0; i < datasets.size(); i++) {
			Dataset dataset = datasets.get(i);
			p = new ScatterplotPanel(1, 1);
			p.setData(Primitives.toPrimitiveArr(dataset.x), Primitives.toPrimitiveArr(dataset.y));
			range = new Range(minpos, -1, maxpos, dataset.getMaxY());
			range.round();
			range.setMinY(0d);
			range.setMinX(44922000);
			p.setDataRange(range);
			p.setTitle(dataset.name);
			p.setPlotElems(true, false);
			grid.addPanel(p);
		}

		grid.draw(output);

		TextFile tfout = new TextFile(output + ".txt", TextFile.W);
		String header = "SNP";
		for (Dataset d : datasets) {
			header += "\t" + d.name;
		}
		tfout.writeln(header);
		for (int s = 0; s < snps.size(); s++) {
			String ln = snps.get(s);
			for (Dataset d : datasets) {
				ln += "\t" + d.y.get(s);
			}
			tfout.writeln(ln);
		}
		tfout.close();

	}


	public void plotIndividualDatasetsGTExFormat(String eqtlgenformat, String[] gtexds, String[] gtexnames, String output) throws IOException, DocumentException {
		TextFile tf = new TextFile(eqtlgenformat, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		int maxpos = -1;
		int minpos = Integer.MAX_VALUE;
		Dataset meta = new Dataset();
		meta.name = "MetaAnalysis";
		while (elems != null) {
			if (elems.length > 1) {
				double metap = Double.parseDouble(elems[0]);
				String snp = elems[1];
				int snpchr = Integer.parseInt(elems[2]);
				int snppos = Integer.parseInt(elems[3]);
				if (snppos > maxpos) {
					maxpos = snppos;
				}
				if (snppos < minpos) {
					minpos = snppos;
				}
//				double metaz = Double.parseDouble(elems[10]);

				meta.x.add((double) snppos);
				double dsp = -Math.log10(metap);
				meta.y.add(dsp);

			}

			ctr++;
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		// now parse GTEx shizzle
		ArrayList<Dataset> datasets = new ArrayList<Dataset>();

		for (int i = 0; i < gtexds.length; i++) {
			Dataset ds = new Dataset();
			ds.name = gtexnames[i];
			TextFile tfg = new TextFile(gtexds[i], TextFile.R);
			String[] gelems = tfg.readLineElems(TextFile.tab);
			while (gelems != null) {
				String snp = gelems[1];
				String[] snpelems = snp.split("_");
				Double pos = Double.parseDouble(snpelems[1]);
				try {
					Double pval = -Math.log10(Double.parseDouble(gelems[6]));

					ds.x.add(pos);
					ds.y.add(pval);
				} catch (NumberFormatException e) {

				}
				gelems = tfg.readLineElems(TextFile.tab);
			}
			tfg.close();
			datasets.add(ds);
		}

		int nrcol = 2;
		int nrrows = (datasets.size()) / nrcol;
		nrrows++;

		Grid grid = new Grid(500, 300, nrrows, 2, 100, 100);

		ScatterplotPanel p = new ScatterplotPanel(1, 1);
		p.setData(Primitives.toPrimitiveArr(meta.x), Primitives.toPrimitiveArr(meta.y));
		Range range = new Range(minpos, 0, maxpos, meta.getMaxY());
		range.round();
		range.setMinY(0d);
		range.setMinX(44922000);
		p.setDataRange(range);
		p.setTitle(meta.name);
		p.setPlotElems(true, false);
		grid.addPanel(p, 0, 0);
		SpacerPanel spacer = new SpacerPanel(1, 1);
		for (int r = 1; r < nrrows; r++) {
			grid.addPanel(spacer, r, 0);
		}

		for (int i = 0; i < datasets.size(); i++) {
			Dataset dataset = datasets.get(i);
			p = new ScatterplotPanel(1, 1);
			p.setData(Primitives.toPrimitiveArr(dataset.x), Primitives.toPrimitiveArr(dataset.y));
			range = new Range(minpos, -1, maxpos, dataset.getMaxY());
			range.round();
			range.setMinY(0d);
			range.setMinX(44922000);
			p.setDataRange(range);
			p.setTitle(dataset.name);
			p.setPlotElems(true, false);
			grid.addPanel(p);
		}

		grid.draw(output);


	}
}
