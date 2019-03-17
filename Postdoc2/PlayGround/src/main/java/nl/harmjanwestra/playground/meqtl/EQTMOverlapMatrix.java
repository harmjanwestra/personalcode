package nl.harmjanwestra.playground.meqtl;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import nl.harmjanwestra.playground.legacy.bedfile.BedFileReader;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.features.FeatureComparator;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.*;

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
	
	public void run2(String featurefileloc, String eqtmfile, String geneAnnotationGTF, String outfolder, int cgwindowsize, boolean exonoverlap) throws IOException {
		
		// load list of annotation files
		
		
		// load a list of cg/gene pairs
		GTFAnnotation gtf = new GTFAnnotation(geneAnnotationGTF);
		HashMap<String, Gene> strToGene = gtf.getStrToGene();
		ArrayList<Pair<Feature, Gene>> eqtms = new ArrayList<>();
		
		QTLTextFile t = new QTLTextFile(eqtmfile, TextFile.R);
		Iterator<EQTL> it = t.getEQtlIterator();
		
		HashSet<String> genehash = new HashSet<String>();
		ArrayList<Gene> genelust = new ArrayList<>();
		HashSet<String> cghash = new HashSet<String>();
		ArrayList<Feature> cglist = new ArrayList<>();
		
		while (it.hasNext()) {
			
			EQTL e = it.next();
			Feature cg = new Feature(Chromosome.parseChr("" + e.getRsChr()), e.getRsChrPos(), e.getRsChrPos());
			cg.setName(e.getRsName());
			
			Gene g = strToGene.get(e.getProbe());
			if (g != null) {
				Pair<Feature, Gene> eqtm = new Pair<>(cg, g);
				eqtms.add(eqtm);
				if (!genehash.contains(e.getProbe())) {
					genelust.add(g);
					genehash.add(e.getProbe());
				}
				if (!cghash.contains(e.getRsName())) {
					cglist.add(cg);
					cghash.add(e.getRsName());
				}
			}
			
			
		}
		t.close();
		
		
		TextFile matrixGenesOut = new TextFile(outfolder + "Overlap_Gene.txt", TextFile.W);
		TextFile matrixCGsOut = new TextFile(outfolder + "Overlap_CG.txt", TextFile.W);
		TextFile matrixEQTMOut = new TextFile(outfolder + "Overlap_EQTM.txt", TextFile.W);
		
		String headerg = "Source\tFeature";
		String headerc = "Source\tFeature";
		String headere = "Source\tFeature";
		
		for (int i = 0; i < eqtms.size(); i++) {
			headere += "\t" + eqtms.get(i).getLeft().getName() + "\t" + eqtms.get(i).getRight().getName();
		}
		for (int i = 0; i < genelust.size(); i++) {
			headerg += "\t" + genelust.get(i).getName();
		}
		for (int i = 0; i < cglist.size(); i++) {
			headerc += "\t" + cglist.get(i);
		}
		matrixCGsOut.writeln(headerc);
		matrixEQTMOut.writeln(headere);
		matrixGenesOut.writeln(headerg);
		
		
		// iterate the list of annotations
		TextFile tf = new TextFile(featurefileloc, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		BedFileReader reader = new BedFileReader();
		while (elems != null) {
			ArrayList<Feature> annotation = reader.readAsList(elems[2]);
			TreeSet<Feature> ts = new TreeSet<>(new FeatureComparator(true));
			ts.addAll(annotation);
			
			String lnoute = elems[0] + "\t" + elems[1];
			String lnoutg = elems[0] + "\t" + elems[1];
			String lnoutc = elems[0] + "\t" + elems[1];
			
			// now iterate each eqtm, cg, and gene
			for (int i = 0; i < eqtms.size(); i++) {
				Pair<Feature, Gene> eqtm = eqtms.get(i);
				Feature cg = eqtm.getLeft();
				Gene gene = eqtm.getRight();
				lnoute += "\t" + overlapcg(ts, cg, cgwindowsize) + "\t" + overlapgene(ts, gene, exonoverlap);
			}
			matrixEQTMOut.writeln(lnoute);
			
			for (int i = 0; i < genelust.size(); i++) {
				Gene gene = genelust.get(i);
				lnoutg += "\t" + overlapgene(ts, gene, exonoverlap);
			}
			matrixGenesOut.writeln(lnoutg);
			
			for (int i = 0; i < cglist.size(); i++) {
				Feature cg = cglist.get(i);
				lnoutc += "\t" + overlapcg(ts, cg, cgwindowsize);
			}
			matrixCGsOut.writeln(lnoutc);
			
			
		}
		
		tf.close();
		
		
		matrixGenesOut.close();
		matrixCGsOut.close();
		matrixEQTMOut.close();
		
	}
	
	private int overlapgene(TreeSet<Feature> ts, Gene gene, boolean exonic) {
		
		Feature fsta = new Feature(gene.getChromosome(), gene.getStart(), gene.getStart());
		Feature fsto = new Feature(gene.getChromosome(), gene.getStop(), gene.getStop());
		
		SortedSet<Feature> subset = ts.subSet(fsta, fsto);
		if (!subset.isEmpty()) {
			return 1;
		}
		
		return 0;
	}
	
	private int overlapcg(TreeSet<Feature> ts, Feature cg, int windowsize) {
		Feature fsta = new Feature(cg.getChromosome(), cg.getStart() - windowsize, cg.getStart() - windowsize);
		Feature fsto = new Feature(cg.getChromosome(), cg.getStop() + windowsize, cg.getStop() + windowsize);
		
		SortedSet<Feature> subset = ts.subSet(fsta, fsto);
		if (!subset.isEmpty()) {
			return 1;
		}
		
		return 0;
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
//		DoubleMatrix2D MeanAndVariance = ds.getMatrix();
//		cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra al = new DenseDoubleAlgebra();
//		MeanAndVariance = al.transpose(MeanAndVariance);
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
