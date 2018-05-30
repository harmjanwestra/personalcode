package nl.harmjanwestra.playground.meqtl;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class CompareCorrelations {
	
	public static void main(String[] args) {
		CompareCorrelations c = new CompareCorrelations();
		String correlationdata = "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2018-05-15-RNACpGCorrelationComparison\\correlationData.txt.gz";
		String plotout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2018-05-15-RNACpGCorrelationComparison\\plot.png";
		
		try {
			c.plot(correlationdata, plotout);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
		System.exit(0);
		
		if (args.length < 6) {
			System.out.println("Usage rna rnasamples meth methsamples eqtms outputfolder");
		} else {
			
			
			try {
				c.run(args[0],
						args[1],
						args[2],
						args[3],
						args[4],
						args[5]);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		
	}
	
	public void plot(String correlationData, String plotout) throws IOException, DocumentException {
		ArrayList<Double> x = new ArrayList<>();
		ArrayList<Double> y = new ArrayList<>();
		
		TextFile tf = new TextFile(correlationData, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length >= 4) {
				Double x1 = Double.parseDouble(elems[1]);
				Double y1 = Double.parseDouble(elems[3]);
				x.add(x1);
				y.add(y1);
//				if (x.size() > 10000) {
//					break;
//				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println(x.size() + " data points loaded");
		
		Grid g = new Grid(500, 500, 1, 1, 100, 100);
		ScatterplotPanel p = new ScatterplotPanel(1, 1);
		p.setAlpha(0.1f);
		p.setData(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
		p.setDataRange(new Range(-1, -1, 1, 1));
		p.setLabels("RNA-correl", "CpG-correl");
		p.setPlotElems(true, false);
		g.addPanel(p);
		
		g.draw(plotout);
		
		
	}
	
	public void run(String rna,
					String rnaSampleLimit,
					String meth,
					String methSampleLimit,
					String eqtmfile,
					String outputfolder) throws Exception {
		
		// load eqtms
		HashSet<String> genes = new HashSet<String>();
		HashSet<String> cpgs = new HashSet<String>();
		
		QTLTextFile tf = new QTLTextFile(eqtmfile, TextFile.R);
		ArrayList<EQTL> eqtms = tf.readList();
		tf.close();
		HashMap<String, String> eqtmlinks = new HashMap<>();
		
		HashMap<String, Boolean> flipCpG = new HashMap<String, Boolean>();
		for (EQTL e : eqtms) {
			String gene = e.getProbe();
			String cpg = e.getRsName();
			if (!genes.contains(gene) && !cpgs.contains(cpg)) {
				if (e.getZscore() < 0) {
					flipCpG.put(cpg, true);
				} else {
					flipCpG.put(cpg, false);
				}
				eqtmlinks.put(gene, cpg);
				genes.add(e.getProbe());
				cpgs.add(e.getRsName());
			}
		}
		
		System.out.println("EQTM file:");
		System.out.println(genes.size() + " genes loaded.");
		System.out.println(cpgs.size() + " cpgs loaded.");
		
		DoubleMatrixDataset<String, String> rnacormat = createCorrelationMatrix(rna, rnaSampleLimit, genes, null);
		
		
		DoubleMatrixDataset<String, String> methcormat = createCorrelationMatrix(meth, methSampleLimit, cpgs, flipCpG);
		
		int[] geneTpCpgIndex = new int[rnacormat.rows()];
		int overlap = 0;
		for (int r = 0; r < rnacormat.rows(); r++) {
			String gene = rnacormat.getRowObjects().get(r);
			String matchingCg = eqtmlinks.get(gene);
			if (matchingCg != null) {
				Integer indexInCorMat = methcormat.getHashRows().get(matchingCg);
				if (indexInCorMat != null) {
					geneTpCpgIndex[r] = indexInCorMat;
					
					overlap++;
				} else {
					geneTpCpgIndex[r] = -1;
				}
			}
		}
		
		if (overlap == 0) {
			System.out.println("ERROR: no overlap between cpgs and genes after correlation matrix");
			System.exit(-1);
		}
		
		TextFile output = new TextFile(outputfolder + "correlationData.txt.gz", TextFile.W);
		// plot a plot maybe?
		Grid g = new Grid(500, 500, 1, 1, 100, 100);
		ScatterplotPanel p = new ScatterplotPanel(1, 1);
		ArrayList<Double> x = new ArrayList<>();
		ArrayList<Double> y = new ArrayList<>();
		
		for (int i = 0; i < rnacormat.rows(); i++) {
			String gene = rnacormat.getRowObjects().get(i);
			String cpg = eqtmlinks.get(gene);
			int icpg = geneTpCpgIndex[i];
			if (icpg >= 0) {
				for (int j = i + 1; j < rnacormat.rows(); j++) {
					int jcpg = geneTpCpgIndex[j];
					if (icpg > -1 && jcpg > -1) {
						// produce output
						String gene2 = rnacormat.getRowObjects().get(j);
						double corrna = rnacormat.getMatrix().getQuick(i, j);
						double corcpg = methcormat.getMatrix().getQuick(icpg, jcpg);
						x.add(corrna);
						y.add(corcpg);
						output.writeln(gene + "\t" + cpg + "\t" + gene2 + "\t" + corrna + "\t" + corcpg);
					}
				}
			}
		}
		
		System.out.println("Now drawing: " + x.size() + " dots");
		p.setAlpha(0.1f);
		p.setData(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
		p.setDataRange(new Range(-1, -1, 1, 1));
		p.setLabels("RNA-correl", "CpG-correl");
		g.addPanel(p);
		g.draw(outputfolder + "outputfile.png");
		
		output.close();
		
	}
	
	private DoubleMatrixDataset<String, String> createCorrelationMatrix(String datasetFile, String colLimit, HashSet<String> rowLimit, HashMap<String, Boolean> rowFlip) throws Exception {
		System.out.println("Calculating correlations for: " + datasetFile);
		HashSet<String> cols = new HashSet<String>();
		TextFile tf = new TextFile(colLimit, TextFile.R);
		cols.addAll(tf.readAsArrayList());
		tf.close();
		
		System.out.println(cols.size() + " samples to select");
		System.out.println(rowLimit.size() + " rows to select");
		
		// load data
		DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadSubsetOfTextDoubleData(datasetFile, "\t", rowLimit, cols);
		System.out.println("Final size: " + ds.rows() + " x " + ds.columns());
		
		DoubleMatrix2D mat = ds.getMatrix();
		
		System.out.println("Z-transform");
		// Z-transform per gene
		for (int r = 0; r < ds.rows(); r++) {
			DoubleMatrix1D row = mat.viewRow(r);
			double[] rowdata = row.toArray();
			String rowStr = ds.getRowObjects().get(r);
			
			boolean flip = false;
			if (rowFlip != null) {
				if (rowFlip.containsKey(rowStr)) {
					flip = rowFlip.get(rowStr);
				}
			}
			
			double mu = Descriptives.mean(rowdata);
			double si = Math.sqrt(Descriptives.variance(rowdata));
			for (int c = 0; c < rowdata.length; c++) {
				double z = (rowdata[c] - mu) / si;
				if (flip) {
					mat.setQuick(r, c, -z);
				} else {
					mat.setQuick(r, c, z);
				}
			}
		}
		
		
		DenseDoubleMatrix2D corrmat = new DenseDoubleMatrix2D(ds.rows(), ds.rows()); // is there a triangular form of this matrix class? should be easy to represent in such a format
		System.out.println("Output is " + ds.rows() + " x " + ds.rows());
		ProgressBar pb = new ProgressBar(ds.rows(), "Correlating!");
		for (int i = 0; i < ds.rows(); i++) {
			corrmat.setQuick(i, i, 1);
			DoubleMatrix1D rowi = mat.viewRow(i);
			double[] rowdatai = rowi.toArray();
			for (int j = i + 1; j < ds.rows(); j++) {
				DoubleMatrix1D rowj = mat.viewRow(j);
				double[] rowdataj = rowj.toArray();
				
				double r = JSci.maths.ArrayMath.correlation(rowdatai, rowdataj);
				corrmat.setQuick(i, j, r);
				corrmat.setQuick(j, i, r);
			}
			pb.iterate();
		}
		pb.close();
		
		System.out.println("Done. ");
		DoubleMatrixDataset<String, String> output = new DoubleMatrixDataset<>();
		output.setMatrix(corrmat);
		output.setRowObjects(ds.getRowObjects());
		output.setColObjects(ds.getRowObjects());
		
		return output;
	}
}
