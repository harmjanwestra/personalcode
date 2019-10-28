package nl.harmjanwestra.playground.biogen;

import com.itextpdf.text.DocumentException;
import umcg.genetica.containers.Pair;
import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.panels.Panel;
import umcg.genetica.graphics.themes.DefaultTheme;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

import java.awt.*;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Objects;

public class ForestPlot {

	public static void main(String[] args) {
		String file = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Cis-Genes\\Primary\\eQTLsFDR0.05-ProbeLevel.txt.gz";

		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\GSTO2\\ForestPlot.pdf";

		String gene = "ENSG00000065621.14";
		String snp = "10:104279545:rs276203:A_G";

		ForestPlot p = new ForestPlot();
		try {
			p.run(file, output, gene, snp);

		} catch (IOException | DocumentException e) {
			e.printStackTrace();
		}
	}

	public void run(String efile, String output, String qgene, String qsnp) throws IOException, DocumentException {
		TextFile tf = new TextFile(efile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			String snp = elems[1];
			String gene = elems[4];

			if ((qsnp == null && qgene == null) || (qgene.equals(gene) && qsnp.equals(snp))) {
				System.out.println("Found eQTL, pval "+elems[0]);

				String betastr = elems[19];
				String nstr = elems[13];
				String namestr = elems[11];
				String metabstr = elems[18];

				String[] betaelems = betastr.split(";");
				String[] nelems = nstr.split(";");
				String[] namesarr = namestr.split(";");

				ArrayList<Double> betas = new ArrayList<>();
				ArrayList<Double> ses = new ArrayList<>();
				ArrayList<Integer> ns = new ArrayList<>();
				ArrayList<String> names = new ArrayList<>();

				int sum = 0;
				for (int i = 0; i < betaelems.length; i++) {
					Pair<Double, Double> q = parsebeta(betaelems[i]);
					betas.add(-q.getLeft());
					ses.add(q.getRight());
					try {
						int n = Integer.parseInt(nelems[i]);
						ns.add(n);
						sum += n;

					} catch (NumberFormatException e) {
						ns.add(0);
					}
					;
					names.add(namesarr[i]);

				}

//		ns.add(sum);

				Pair<Double, Double> metabetavals = parsebeta(metabstr);

//		names.add("MetaBrain");


				ForestplotPanel p = new ForestplotPanel(1, 1);
				p.setData(Primitives.toPrimitiveArr(betas),
						Primitives.toPrimitiveArr(ses),
						Primitives.toPrimitiveArr(ns),
						names.toArray(new String[0]),
						-metabetavals.getLeft(),
						metabetavals.getRight(),
						sum


				);

				Grid g = new Grid(250, 750, 1, 1, 300, 300);
				g.addPanel(p);
				g.draw(output + gene + ".pdf");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

	}

	private Pair<Double, Double> parsebeta(String betaelem) {

		String[] betaelems = betaelem.split(" ");
		if (betaelems.length > 1) {
			Double beta = Double.parseDouble(betaelems[0]);
			String sestr = betaelems[1].replaceAll("\\(", "");
			sestr = sestr.replaceAll("\\)", "");
			Double se = Double.parseDouble(sestr);
			return new Pair<Double, Double>(beta, se);


		}
		return new Pair<Double, Double>(Double.NaN, Double.NaN);
	}

	public class ForestplotPanel extends Panel {

		private double metabeta;
		private double metase;
		private int metan;
		private double maxeffectsize = 1d;


		public ForestplotPanel(int nrRows, int nrCols) {
			super(nrRows, nrCols);
		}

		public void setMaxeffectsize(double beta) {
			this.maxeffectsize = maxeffectsize;
		}

		double[] beta = null;
		double[] se = null;
		int[] n = null;
		String[] names = null;

		public void setData(double[] beta, double[] se, int[] n, String[] names, double metabeta, double metase, int metan) {
			this.beta = beta;
			this.se = se;
			this.n = n;
			this.names = names;

			this.metabeta = metabeta;
			this.metase = metase;
			this.metan = metan;

		}

		private class Dataset {
			String name;
			double beta;
			double se;
			int n;

			@Override
			public boolean equals(Object o) {
				if (this == o) return true;
				if (o == null || getClass() != o.getClass()) return false;
				Dataset dataset = (Dataset) o;
				return Double.compare(dataset.beta, beta) == 0;
			}

			@Override
			public int hashCode() {
				return Objects.hash(beta);
			}
		}

		private class DatasetComparator implements Comparator<Dataset> {

			@Override
			public int compare(Dataset o1, Dataset o2) {
				if (o1.equals(o2)) {
					return 0;
				} else if (o1.beta > o2.beta) {
					return -1;
				} else {
					return 1;
				}

			}
		}

		@Override
		public void draw(DefaultGraphics defaultGraphics) {


			int maxn = 0;
			int minn = 100000;
			ArrayList<Dataset> datasets = new ArrayList<>();
			for (int i = 0; i < n.length; i++) {
				if (n[i] > maxn) {
					maxn = n[i];
				}
				if (n[i] < minn) {
					minn = n[i];
				}
				Dataset s = new Dataset();
				s.beta = beta[i];
				s.n = n[i];
				s.name = names[i];
				s.se = se[i];
				datasets.add(s);
			}

			Collections.sort(datasets, new DatasetComparator());

			int maxsize = 50;

			int halfwidht = width / 2;

			Graphics2D g2d = defaultGraphics.getG2d();

			DefaultTheme theme = new DefaultTheme();
			g2d.setColor(theme.getLightGrey());

			g2d.setStroke(theme.stroke);


			// plot range
			g2d.drawString("" + maxeffectsize, x0 + width, y0 - 10);
			g2d.drawString("" + (-maxeffectsize), x0, y0 - 10);
			g2d.drawString("0", x0 + (width / 2), y0 - 10);

			// plot header
			g2d.drawString("Dataset", x0 + width + 20, y0 - 10);
			g2d.drawString("N", x0 + width + 150, y0 - 10);
			g2d.drawString("Beta", x0 + width + 200, y0 - 10);
			g2d.drawString("Pval", x0 + width + 250, y0 - 10);

			// ------
			g2d.drawLine(x0, y0, x0 + width, y0);
			g2d.setStroke(theme.strokeDashed);
			//    |
			g2d.drawLine(x0 + halfwidht, y0, x0 + halfwidht, y0 + height);
			g2d.setStroke(theme.stroke);
			// ------
			g2d.drawLine(x0, y0 + height, x0 + width, y0 + height);


			int marginbetweenbox = 40;
			int boxsizemax = 25;
			int boxsizemin = 10;

			int starty = y0 + marginbetweenbox;

			for (int i = 0; i < datasets.size(); i++) {

				Dataset d = datasets.get(i);

				plotbox(d, boxsizemin, boxsizemax, minn, maxn, halfwidht, maxeffectsize, g2d, starty);

				starty += marginbetweenbox;
			}

			Dataset metad = new Dataset();
			metad.beta = metabeta;
			metad.se = metase;
			metad.n = metan;
			metad.name = "MetaBrain";
			plotbox(metad, boxsizemin, boxsizemax, minn, maxn, halfwidht, maxeffectsize, g2d, starty);


		}

		private void plotbox(Dataset d, int boxsizemin, int boxsizemax, int minn, int maxn, int halfwidht, double range, Graphics2D g2d, int starty) {
			double b = d.beta;
			if (!Double.isNaN(b)) {

				int boxsize = (int) Math.floor(boxsizemin + ((double) (d.n - minn) / (maxn - minn)) * (boxsizemax - boxsizemin));

				int x = x0 + halfwidht;
				if (b > 0) {
					if (b > range) {
						b = range;
					}
					int nrpx = (int) Math.floor((b / range) * halfwidht);
					x += nrpx;
				} else if (b < 0) {
					if (b < -range) {
						b = -range;
					}
					int nrpx = (int) Math.floor((b / range) * halfwidht);
					x += nrpx;
				}

				g2d.fillRect(x - (boxsize / 2), starty, boxsize, boxsize);
				int nrpxse = (int) Math.floor((d.se / range) * halfwidht);
				int minpx = x - nrpxse;
				int maxpx = x + nrpxse;
				g2d.drawLine(minpx, starty + (boxsize / 2), maxpx, starty + (boxsize / 2));

				g2d.drawString(d.name, x0 + width + 20, starty);
				g2d.drawString("" + d.n, x0 + width + 150, starty);
				double z = (d.beta / d.se);
				double pval = ZScores.zToP(z);
				DecimalFormat df = new DecimalFormat("0.0E000");
				DecimalFormat df2 = new DecimalFormat("0.000");
				String pvalstr = df.format(pval);
				if (pval >= 0.0001) {
					pvalstr = df2.format(pval);
				}

				System.out.println(d.name + "\t" + pvalstr + "\t" + pval);


				g2d.drawString("" + df2.format(d.beta), x0 + width + 200, starty);
				g2d.drawString("" + pvalstr, x0 + width + 250, starty);

			}
		}
	}
}
