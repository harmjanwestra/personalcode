package nl.harmjanwestra.playground.biogen.locusdebug;

import com.itextpdf.text.DocumentException;

import umcg.genetica.containers.Triple;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.graphics.panels.SpacerPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Objects;

public class LinRegOutlierIdentify {


	public static void main(String[] argS) {
		String input = "C:\\Users\\harm-jan\\Downloads\\zsore_table.topSNPs.txt";
		LinRegOutlierIdentify l = new LinRegOutlierIdentify();
		String outdir = "C:\\Users\\harm-jan\\Downloads\\tmpout\\";

		try {
			l.run(input, outdir);
		} catch (IOException | DocumentException e) {
			e.printStackTrace();
		}
	}


	public void run(String input, String outdir) throws IOException, DocumentException {


		ArrayList<Dataset> datasets = new ArrayList<>();
		TextFile tf = new TextFile(input, TextFile.R);

		String[] header = tf.readLineElems(TextFile.tab);
		for (int i = 5; i < header.length; i++) {
			Dataset d = new Dataset();
			d.name = header[i];
			datasets.add(d);
		}

		// load
		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<String> eqtls = new ArrayList<>();
		int ctr = 0;
		while (elems != null) {

			String eqtl = elems[0] + ";" + elems[1];
			eqtls.add(eqtl);
			for (int i = 5; i < elems.length; i++) {
				Double v = null;
				try {
					v = Double.parseDouble(elems[i]);
				} catch (NumberFormatException e) {

				}
				datasets.get(i - 5).vals.add(v);
			}
			ctr++;
			if (ctr % 1000 == 0) {
				System.out.println(ctr + " read.");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		// compare datasets

		Grid g = new Grid(500, 500, datasets.size(), datasets.size(), 100, 100);
		for (int i = 0; i < datasets.size(); i++) {
			Dataset d1 = datasets.get(i);
			for (int j = 0; j < i + 1; j++) {
				g.addPanel(new SpacerPanel(1, 1));
			}
			for (int j = i + 1; j < datasets.size(); j++) {
				Dataset d2 = datasets.get(j);
				Triple<ArrayList<Double>, ArrayList<Double>, ArrayList<String>> intersect = d1.intersect(d2, eqtls);

				// sort by d1?
				// very dumb, but wrap eqtls in some kind of eqtl class

				ArrayList<Double> di1 = intersect.getLeft();
				ArrayList<Double> di2 = intersect.getMiddle();
				ArrayList<String> ei = intersect.getRight();

				ArrayList<eqtl> compeqtls = new ArrayList<>();
				for (int q = 0; q < di1.size(); q++) {
					eqtl e = new eqtl();
					e.name = ei.get(q);
					e.v1 = di1.get(q);
					e.v2 = di2.get(q);
					compeqtls.add(e);
				}

//				System.out.println("Test !");
//				eqtlcomparator<eqtl> comp = new eqtlcomparator<eqtl>();
//				for (int q = 0; q < compeqtls.size(); q++) {
//					for (int r = q + 1; r < compeqtls.size(); r++) {
//						int c1 = comp.compare(compeqtls.get(q), compeqtls.get(r));
//						int c2 = comp.compare(compeqtls.get(r), compeqtls.get(q));
//						boolean equal1 = compeqtls.get(q).equals(compeqtls.get(r));
//						boolean equal2 = compeqtls.get(r).equals(compeqtls.get(q));
//
//
//						if ((c1 == 0 && c2 != 0) || (c2 == 0 && c1 != 0)) {
//							System.out.println("Error 1: " + compeqtls.get(q) + "\t" + compeqtls.get(r) + "\t" + c1 + "\t" + c2 + "\t" + equal1 + "\t" + equal2);
//						} else if ((c1 < 0 && c2 > 0) || (c1 > 0 && c2 < 0)) {
//							System.out.println("Error 2: " + compeqtls.get(q) + "\t" + compeqtls.get(r) + "\t" + c1 + "\t" + c2 + "\t" + equal1 + "\t" + equal2);
//						}
//
//					}
//				}
				System.out.println(i + " vs " + j);
				Collections.sort(compeqtls, new eqtlcomparator<eqtl>());

				int window = 100;
				int stepsize = 1;
				int pos = 0;
				double zthresh = 4;

				while (pos < compeqtls.size()) {
					if (pos + window > compeqtls.size()) {
						window = compeqtls.size() - pos;
					}

					double[] x = new double[window];
					double[] y = new double[window];
					for (int q = 0; q < window; q++) {
						eqtl e = compeqtls.get(q + pos);
						x[q] = e.v1;
						y[q] = e.v2;
					}

					double meanx = Descriptives.mean(x);
					double meany = Descriptives.mean(y);
					double sdx = Math.sqrt(Descriptives.variance(x));
					double sdy = Math.sqrt(Descriptives.variance(y));

					// determine outliers at z>threshold
					for (int q = pos; q < pos + window; q++) {
						eqtl e = compeqtls.get(q);
						double z1 = (e.v1 - meanx) / sdx;
						double z2 = (e.v2 - meany) / sdy;

						if (Math.abs(z1) > zthresh) {
							// outlier
							e.outlierx = true;
						}
						if (Math.abs(z2) > zthresh) {
							// outlier
							e.outliery = true;
						}

					}
					pos += stepsize;
				}


				ArrayList<Double> x1 = new ArrayList<>();
				ArrayList<Double> xoutlier = new ArrayList<>();
				ArrayList<Double> y1 = new ArrayList<>();
				ArrayList<Double> youtlier = new ArrayList<>();

				int nroutliers = 0;
				for (int q = 0; q < compeqtls.size(); q++) {
					eqtl e = compeqtls.get(q);
					if (e.outlierx || e.outliery) {
						xoutlier.add(e.v1);
						youtlier.add(e.v2);
						nroutliers++;
					} else {
						x1.add(e.v1);
						y1.add(e.v2);
					}
				}

				double[][] x = new double[2][];
				double[][] y = new double[2][];
				x[0] = Primitives.toPrimitiveArr(x1);
				y[0] = Primitives.toPrimitiveArr(y1);
				x[1] = Primitives.toPrimitiveArr(xoutlier);
				y[1] = Primitives.toPrimitiveArr(youtlier);


				ScatterplotPanel p1 = new ScatterplotPanel(1, 1);

				p1.setData(x, y);
				p1.setLabels(d1.name, d2.name);
				p1.setPlotElems(true, true);
				p1.setDatasetLabels(new String[]{"Inliers", "Outliers"});
				p1.setTitle(nroutliers + " / " + compeqtls.size() + " outliers.");
				g.addPanel(p1);

//				System.exit(-1);
			}
		}
		g.draw(outdir + "comp.pdf");

	}

	private class eqtl {
		String name;
		double v1;
		double v2;
		boolean outlierx = false;
		boolean outliery = false;

		public boolean equals(eqtl o) {
			return Double.compare(o.v1, v1) == 0;
		}

		@Override
		public int hashCode() {
			return Objects.hash(v1);
		}

		@Override
		public String toString() {
			return "eqtl{" +
					"name='" + name + '\'' +
					", v1=" + v1 +
					", v2=" + v2 +
					'}';
		}
	}


	private class eqtlcomparator<E> implements Comparator<eqtl> {
		@Override
		public int compare(eqtl o1, eqtl o2) {
			if (o1.equals(o2)) {
				return 0;
			}
			if (o1.v1 > o2.v1) {
				return 1;
			} else if (o1.v1 < o2.v1) {
				return -1;
			} else {
				return 0;
			}
		}
	}

	private class Dataset {
		ArrayList<Double> vals = new ArrayList<>();
		String name;

		public Triple<ArrayList<Double>, ArrayList<Double>, ArrayList<String>> intersect(Dataset other, ArrayList<String> eqtls) {
			ArrayList<Double> a1 = new ArrayList<>();
			ArrayList<Double> a2 = new ArrayList<>();
			ArrayList<String> eqtlout = new ArrayList<>();
			for (int i = 0; i < vals.size(); i++) {
				if (vals.get(i) != null && other.vals.get(i) != null) {
					a1.add(vals.get(i));
					a2.add(other.vals.get(i));
					eqtlout.add(eqtls.get(i));
				}
			}
			return new Triple<>(a1, a2, eqtlout);
		}
	}
}
