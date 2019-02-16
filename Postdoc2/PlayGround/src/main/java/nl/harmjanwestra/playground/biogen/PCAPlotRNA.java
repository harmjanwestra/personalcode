package nl.harmjanwestra.playground.biogen;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class PCAPlotRNA {
	
	public static void main(String[] args) {
		
		String[] pcafiles = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-12-03-RNASeqOutliers\\Braineac-pc1_4.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-12-03-RNASeqOutliers\\CMC-pc1_4.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-12-03-RNASeqOutliers\\GTEx-pc1_4.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-12-03-RNASeqOutliers\\MayoCER-pc1_4.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-12-03-RNASeqOutliers\\MayoTCX-pc1_4.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-12-03-RNASeqOutliers\\MSBB-pc1_4.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-12-03-RNASeqOutliers\\ROSMAP-pc1_4.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-12-03-RNASeqOutliers\\TargetALS-pc1_4.txt"
		};
		
		String[] datasetnames = new String[]{"Braineac", "CMC", "GTEx", "MayoCER", "MayoTCX", "MSBB", "ROSMAP", "TargetALS"};
		
		int c1 = 1;
		int c2 = 2;
		String outlierfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-12-03-RNASeqOutliers\\rna-pcaoutliers.txt";
		
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-12-03-RNASeqOutliers\\rna-pcaoutliers.pdf";
		
		PCAPlotRNA p = new PCAPlotRNA();
		try {
			p.run(pcafiles, datasetnames, c1, c2, outlierfile, out);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String[] pcafiles, String[] datasetnames, int c1, int c2, String outlierfile, String out) throws IOException, DocumentException {
		
		Grid grid = new Grid(500, 500, 2, 4, 100, 100);
		
		for (int d = 0; d < pcafiles.length; d++) {
			String pcafile = pcafiles[d];
			String datasetname = datasetnames[d];
			
			HashSet<String> outlierIDs = new HashSet<String>();
			TextFile tf = new TextFile(outlierfile, TextFile.R);
			String ln = tf.readLine();
			while (ln != null) {
				outlierIDs.add(ln);
				ln = tf.readLine();
			}
			tf.close();
			
			
			ArrayList<Pair<Double, Double>> outliers = new ArrayList<Pair<Double, Double>>();
			ArrayList<Pair<Double, Double>> inliers = new ArrayList<Pair<Double, Double>>();
			TextFile tf2 = new TextFile(pcafile, TextFile.R);
			String[] header = tf2.readLineElems(TextFile.tab);
			String[] labels = new String[]{header[c1], header[c2]};
			String[] elems = tf2.readLineElems(TextFile.tab);
			while (elems != null) {
				String sample = elems[0];
				Double x = Double.parseDouble(elems[c1]);
				Double y = Double.parseDouble(elems[c2]);
				Pair<Double, Double> v = new Pair<>(x, y);
				if (outlierIDs.contains(sample)) {
					outliers.add(v);
				} else {
					inliers.add(v);
				}
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();
			
			double[][] x = new double[2][];
			double[][] y = new double[2][];
			x[0] = new double[inliers.size()];
			x[1] = new double[outliers.size()];
			y[0] = new double[inliers.size()];
			y[1] = new double[outliers.size()];
			for (int p = 0; p < inliers.size(); p++) {
				x[0][p] = inliers.get(p).getLeft();
				y[0][p] = inliers.get(p).getRight();
			}
			for (int p = 0; p < outliers.size(); p++) {
				x[1][p] = outliers.get(p).getLeft();
				y[1][p] = outliers.get(p).getRight();
			}
			
			DefaultTheme def = new DefaultTheme();
			Color[] colors = new Color[]{
					new Color(0, 0, 0),
					
					new Color(255, 51, 51),
			};
			def.setColors(colors);
			
			int nroutliers = outliers.size();
			int nrnonoutliers = inliers.size();
			int total = nrnonoutliers + nroutliers;
			
			
			ScatterplotPanel p = new ScatterplotPanel(1, 1);
			p.setTitle(datasetname + " - " + nroutliers + " outliers / " + nrnonoutliers + " remain");
			p.setData(x, y);
			p.setDatasetLabels(new String[]{"Non-outliers", "Outliers"});
			p.setLabels(labels[0], labels[1]);
			p.setPlotElems(true, true);
			p.setTheme(def);
			p.setAlpha(0.8f);
			grid.addPanel(p);
		}
		
		grid.draw(out);
		
		
	}
}
