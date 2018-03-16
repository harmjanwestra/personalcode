package nl.harmjanwestra.playground.meqtl;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

public class EQTMOverlapMatrix {
	
	
	public static void main(String[] args) {
		
		
		String annotdir = "D:\\Work\\Data\\Enhancers\\Roadmap\\consolidatedImputedGappedPeak\\";
		String cgfile = "D:\\Work\\Projects\\2017-11-eQTLMeta\\eqtm\\eQTLsFDR0.0-SNPLevel-flipped.txt.gzsplit.txt";
		String markfilter = "H3K27ac";
		String celltypefilter = "D:\\Work\\Projects\\2018-eQTMPredict\\annotation\\roadmap-celltypes-blood-MJ.txt";
		EQTMOverlapMatrix om = new EQTMOverlapMatrix();
		String outfile = cgfile + markfilter;
		
		try {
			om.run(annotdir, markfilter, celltypefilter, cgfile, outfile);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void run(String annotfileloc, String markfilter, String celltypefilter, String cgfile, String outf) throws IOException {


//		TextFile tf = new TextFile(annotfileloc, TextFile.R);
//		ArrayList<String> fileslist = tf.readAsArrayList();
//		tf.close();
//		System.out.println(fileslist.size() + " files");
		
		TextFile tf = new TextFile(celltypefilter, TextFile.R);
		ArrayList<String> set = tf.readAsArrayList();
		tf.close();
		
		File folder = new File(annotfileloc);
		File[] files = folder.listFiles();
		
		ArrayList<File> selectedFiles = new ArrayList<>();
		for (File f : files) {
			boolean include = false;
			if (f.getName().contains(markfilter)) {
				for (String s : set) {
					if (f.getName().contains(s)) {
						include = true;
					}
				}
			}
			if (include) {
				System.out.println("Including: " + f.getName());
				selectedFiles.add(f);
			}
		}
		System.out.println(selectedFiles.size() + " files match patterns");
		
		ArrayList<Feature> cgs = new ArrayList<>();
		TextFile tf2 = new TextFile(cgfile, TextFile.R);
		tf2.readLine();
		String[] elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			Feature f = new Feature();
			f.setName(elems[1]);
			f.setChromosome(Chromosome.parseChr(elems[2]));
			f.setStart(Integer.parseInt(elems[3]) - 25);
			f.setStop(Integer.parseInt(elems[3]) + 25);
			Double z = Double.parseDouble(elems[10]);
			if (z > 0) {
				f.setNrN(1);
			} else {
				f.setNrN(0);
			}
			cgs.add(f);
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		
		System.out.println(cgs.size() + " cgs");
		
		Collections.sort(cgs, new FeatureComparator(true));
		
		// iterate over the files
		double[] y = new double[cgs.size()];
		for (int i = 0; i < cgs.size(); i++) {
			y[i] = cgs.get(i).getNrN();
		}
		
		
		double[][] matrix = new double[selectedFiles.size()][cgs.size()];
		BedFileReader r = new BedFileReader();
		TextFile out = new TextFile(outf + "-matrix.txt.gz", TextFile.W);
		TextFile outc = new TextFile(outf + "-correl.txt.gz", TextFile.W);
		// concat cg names
		String header = "-";
		for (int c = 0; c < cgs.size(); c++) {
			header += "\t" + cgs.get(c).getName();
		}
		out.writeln(header);
		String header2 = "PosOrNegEQTMCorrelation";
		for (int c = 0; c < cgs.size(); c++) {
			header2 += "\t" + cgs.get(c).getNrN();
		}
		
		for (int q = 0; q < selectedFiles.size(); q++) {
			System.out.println(q + "/" + selectedFiles.size());
			ArrayList<Feature> feats = r.readAsList(selectedFiles.get(q).toPath().toAbsolutePath().toString(), cgs);
			System.out.println(feats.size() + " elements ");
			double[] x = new double[cgs.size()];
//			for (int i = 0; i < feats.size(); i++) {
			for (int j = 0; j < cgs.size(); j++) {
				if (cgs.get(j).overlaps(feats)) {
					matrix[q][j] = 1;
					x[j] = 1;
//					if (cgs.get(j).getNrN() > 0) {
//						nrpos++;
//					} else {
//						nrneg++;
//					}
//
				}
			}
//			}
			PearsonsCorrelation pc = new PearsonsCorrelation();
			double correl = pc.correlation(x, y);
			String outcln = selectedFiles.get(q).getName() + "\t" + correl;
			outc.writeln(outcln);
			
			File fileObj = selectedFiles.get(q);
			String name = fileObj.getName();
			System.out.println(name + "\t" + correl);
			String outln = name + "\t" + Strings.concat(matrix[q], Strings.tab);
			out.writeln(outln);
		}
		out.close();
		outc.close();

//		DoubleMatrixDataset ds = new DoubleMatrixDataset<>();
//		ds.setMatrix(matrix);
//		DoubleMatrix2D mat = ds.getMatrix();
//		cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra al = new DenseDoubleAlgebra();
//		mat = al.transpose(mat);
//
//		File fileObj = new File(fileslist.get(q));
//		String name = fileObj.getName();
//		for (int c = 0; c < cgs.size(); c++) {
//			header += "\t" + cgs.get(c).getName();
//		}
//
//
//		TextFile outd = new TextFile(outf + "-nrposneg.txt", TextFile.W);
//		for (int j = 0; j < cgs.size(); j++) {
//
//		}
//		outd.close();
	
	}
	
}
