package nl.harmjanwestra.playground.biogen.covariates;

import com.itextpdf.text.DocumentException;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.HeatmapPanel;
import umcg.genetica.graphics.panels.HistogramPanel;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.graphics.themes.DefaultTheme;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Descriptives;

import java.awt.*;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Objects;
import java.util.stream.IntStream;

public class CovariatePlotter {

	public static void main(String[] args) {
		String cov = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates-qualityscores-filter.txt";
//        String ds = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna\\Freeze1PlusENA-SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz";
		String ds = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna\\top9variance\\Freeze1PlusENA-SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.txt.gz";

		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\covariatecorrelationwithexp\\correlAfterTop9";


		CovariatePlotter cv = new CovariatePlotter();
		try {
			ds = "D:\\Work\\Freeze2\\run2\\2019-04-11-Freeze2.TMM.SampleFilter.SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed";
			cov = "D:\\Work\\Freeze2\\run2\\2019-04-11-Freeze2.TMM.Covariates-Numeric.txt";
			output = "D:\\Work\\Freeze2\\run2\\2019-04-09-CovCorrelBeforeCorrection-pearson.txt";
			String individualplotprefix = "D:\\Work\\Freeze2\\run2\\plots\\2019-04-11-r1-";

			boolean spearman = false;
			cv.correlateCovariatesWithExp(ds, cov, output, individualplotprefix, spearman, 1.1);


		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void correlateCovariatesWithExp(String inputxyfile, String covariatefile, String outplot, String individualplotprefix, boolean spearman, double corthreshold) throws Exception, DocumentException {

		DefaultTheme def = new DefaultTheme();
		Color[] colors = new Color[]{
				Color.decode("#f44336".toUpperCase()),
				Color.decode("#2196f3".toUpperCase()),
				Color.decode("#ffeb3b".toUpperCase()),
				Color.decode("#e91e63".toUpperCase()),
				Color.decode("#009688".toUpperCase()),
				Color.decode("#8bc34a".toUpperCase()),
				Color.decode("#ff5722".toUpperCase()),

				Color.decode("#3f51b5".toUpperCase()),
				Color.decode("#00bcd4".toUpperCase()),
				Color.decode("#9c27b0".toUpperCase()),
				Color.decode("#673ab7".toUpperCase()),

				Color.decode("#03a9f4".toUpperCase()),
				Color.decode("#607d8b".toUpperCase()),
				Color.decode("#4caf50".toUpperCase()),
				Color.decode("#cddc39".toUpperCase()),
				Color.decode("#ffc107".toUpperCase()),
				Color.decode("#ff9800".toUpperCase()),
				Color.decode("#795548".toUpperCase()),

		};

		int alpha = 80;
		for (int c = 0; c < colors.length; c++) {
			colors[c] = new Color(colors[c].getRed(), colors[c].getGreen(), colors[c].getBlue(), alpha);
		}

		def.setColors(colors);

		Color lg = def.getLightGrey();
		def.setLightgrey(new Color(lg.getRed(), lg.getGreen(), lg.getBlue(), alpha));

		DoubleMatrixDataset<String, String> dscovariate = DoubleMatrixDataset.loadDoubleData(covariatefile); // samples on rows

		HashSet<String> samples = new HashSet<String>();
		samples.addAll(dscovariate.getRowObjects());


		double[][] heatmap = new double[dscovariate.columns()][dscovariate.columns()];
		for (int d = 0; d < dscovariate.columns(); d++) {
			double[] x = dscovariate.getMatrix().viewColumn(d).toArray();
			PCAPlot p = new PCAPlot();
			for (int d2 = 0; d2 < dscovariate.columns(); d2++) {
				double[] y = dscovariate.getMatrix().viewColumn(d).toArray();
				heatmap[d][d2] = p.pruneAndCorrelate(x, y);
			}
		}


//        Grid grid2 = new Grid(500, 500, 1, 1, 100, 150);
//        HeatmapPanel hm = new HeatmapPanel(1, 1);
//        hm.setData(heatmap, dscovariate.getColObjects().toArray(new String[0]), dscovariate.getColObjects().toArray(new String[0]));
//        hm.setPlotMode(HeatmapPanel.MODE.FULL);
//        hm.
//        grid2.addPanel(hm);

//        grid2.draw(outplot + "-cov-hm.pdf");
//        System.exit(-1);

		DoubleMatrixDataset<String, String> ds;
		if (inputxyfile.endsWith(".txt.gz")) {
			ds = DoubleMatrixDataset.loadSubsetOfTextDoubleData(inputxyfile, '\t', null, samples); // samples on rows
		} else {
			ds = DoubleMatrixDataset.loadDoubleBinaryData(inputxyfile);
		}
		int[] samplemap = new int[ds.columns()];
		int nrsharedsamples = 0;
		System.out.println("Missing samples...");
		for (int r = 0; r < ds.columns(); r++) {

			String sample = ds.getColObjects().get(r);

			Integer index = dscovariate.getHashRows().get(sample);
			if (index == null) {
				index = -1;
				System.out.println(sample);
			} else {
				nrsharedsamples++;

			}
			samplemap[r] = index;
		}


		System.out.println(nrsharedsamples + " samples shared. ");
		if (nrsharedsamples == 0) {
			System.exit(-1);
		}
		int nrbins = 10;
		int nrcols = 5;
		Grid grid = new Grid(250, 150, (int) Math.ceil(dscovariate.columns() / nrcols) + 4, nrcols, 100, 150);

		double[][] outputmat = new double[ds.rows()][dscovariate.columns()];

		HistogramObj[] data = new HistogramObj[dscovariate.columns()];

		ProgressBar pb = new ProgressBar(dscovariate.columns());
		int finalNrsharedsamples = nrsharedsamples;
		HistogramObj[] finalData = data;
		IntStream.range(0, dscovariate.columns()).parallel().forEach(c -> {
//        IntStream.range(0, 1).parallel().forEach(c -> {
			String cov = dscovariate.getColObjects().get(c);
			int[] bins = new int[nrbins];
			double[] covariate = dscovariate.getCol(c).toArray();

			double[] allcorrel = new double[ds.rows()];
			String covariatename = dscovariate.getColObjects().get(c);
			double[] finalcovariate = new double[ds.columns()];
			int nonmissingvalues = 0;
			for (int i = 0; i < ds.columns(); i++) {
				int map = samplemap[i];
				if (map == -1) {
					finalcovariate[i] = Double.NaN;
				} else {
					double v = covariate[map];
					finalcovariate[i] = v;
					if (!Double.isNaN(v)) {
						nonmissingvalues++;
					}
				}
			}

			if (nonmissingvalues != finalNrsharedsamples) {
				System.out.println("Excluding covar " + cov + " because n=" + nonmissingvalues);
			} else if (nonmissingvalues == finalNrsharedsamples) {
				PCAPlot p = new PCAPlot();

				double maxcor = 0;
				String maxgene = null;
				DecimalFormat format = new DecimalFormat("#.###");
				for (int g = 0; g < ds.rows(); g++) {
					double[] x = ds.getRow(g).toArray();
					String gene = ds.getRowObjects().get(g);
					double xcor = 0;
					try {
						xcor = p.pruneAndCorrelate(x, finalcovariate, corthreshold, individualplotprefix, gene, cov, spearman);
					} catch (IOException e) {
						e.printStackTrace();
					} catch (DocumentException e) {
						e.printStackTrace();
					}
					double rsq = xcor * xcor;


					if (Math.abs(rsq) > maxcor) {
						maxcor = rsq;
						maxgene = ds.getRowObjects().get(g);
					}

					int bin = (int) Math.floor(rsq * nrbins);
					if (bin > nrbins - 1) {
						bin = nrbins - 1;
					}

					allcorrel[g] = rsq;
					bins[bin]++;
					outputmat[g][c] = xcor;
				}

				double meanrsq = Descriptives.mean(allcorrel);
				double varrsq = Descriptives.variance(allcorrel);

				HistogramPanel hist = new HistogramPanel(1, 1);
				hist.setData(bins);
				hist.setAxisLabels("Rsq", "Freq");
				hist.setTitle(covariatename + " max rsq: " + format.format(maxcor));
				hist.setTheme(def);
				hist.setDatasetPlotTypes(new HistogramPanel.DATASETPLOTTYPE[]{HistogramPanel.DATASETPLOTTYPE.BAR});
				hist.setMarginBetweenBinClusters(0);
				HistogramObj o = new HistogramObj();
				hist.setMarginBetweenBinClusters(0);
				hist.setBinLabels(new String[]{"" + 0, "" + 0.1, "" + 0.2, "" + 0.3, "" + 0.4, "" + 0.5, "" + 0.6, "" + 0.7, "" + 0.8, "" + 0.9});
				hist.setAxisLabels("Spearman R-square", "# of genes");
				o.d = maxcor;
				o.p = hist;
				o.name = covariatename;
				o.gene = maxgene;
				o.varrsq = varrsq;
				o.meanrsq = meanrsq;
				finalData[c] = o;


			}
			pb.iterateSynched();
		});
		pb.close();

		ArrayList<HistogramObj> datatmp = new ArrayList<>();
		for (int c = 0; c < finalData.length; c++) {
			if (finalData[c] != null) {
				datatmp.add(data[c]);
			}
		}

		data = datatmp.toArray(new HistogramObj[0]);
		Arrays.sort(data);
		for (int c = 0; c < data.length; c++) {
			if (data[c] != null) {
				grid.addPanel(data[c].p);
			}
		}

		grid.draw(outplot + ".pdf");

		TextFile tf = new TextFile(outplot + "-max.txt", TextFile.W);
		tf.writeln("Covariate\tgene\tMaxRSq\tMeanRSq\tVarRsq\tZRSq");
		for (int c = 0; c < data.length; c++) {
			tf.writeln(data[c].name + "\t" + data[c].gene + "\t" + data[c].d + "\t" + data[c].meanrsq + "\t" + data[c].varrsq + "\t" + ((data[c].d - data[c].meanrsq) / data[c].varrsq));
		}
		tf.close();

		DoubleMatrixDataset<String, String> outmatObj = new DoubleMatrixDataset<String, String>();
		outmatObj.setMatrix(outputmat);
		outmatObj.setColObjects(dscovariate.getColObjects());
		outmatObj.setRowObjects(ds.getRowObjects());
		outmatObj.save(outplot + ".txt");

	}

	public class HistogramObj implements Comparable<HistogramObj> {

		public HistogramPanel p;
		public String gene;
		public double meanrsq;
		public double varrsq;
		double d;
		String name;

		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			HistogramObj that = (HistogramObj) o;
			return Double.compare(that.d, d) == 0;
		}

		@Override
		public int hashCode() {
			return Objects.hash(d);
		}


		@Override
		public int compareTo(HistogramObj o) {
			if (this.equals(o)) {
				return 0;
			} else if (o.d > this.d) {
				return 1;
			} else {
				return -1;
			}
		}
	}
}
