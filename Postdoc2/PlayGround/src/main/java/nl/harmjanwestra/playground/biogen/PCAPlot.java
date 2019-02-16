package nl.harmjanwestra.playground.biogen;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class PCAPlot {
	
	public static void main(String[] args) {
		String knownClasses = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p3v5a-sampleinfo.txt";
		String superclasses = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p3v5a-superpopulations.txt";
		String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\genotypepca\\plink.eigenvec.txt";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\genotypepca\\AMPAD-";
//		String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\genotypepca\\plink.eigenvec";
//		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\genotypepca\\GTEx-";
//        String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\plink.eigenvec.txt";
//        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-";
//        String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\genotypepca\\plink.eigenvec.txt";
//        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\genotypepca\\CMC-";
//        String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\dnapca\\plink.eigenvec.txt";
//        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\dnapca\\Braineac-";
		
		
		PCAPlot p = new PCAPlot();
		
		try {
			p.plot(pcafile, knownClasses, superclasses, output);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
		
	}
	
	
	public void plot(String pcafile, String knownclasses, String superclasses, String output) throws IOException, DocumentException {
		
		HashMap<String, String> superpop = new HashMap<>();
		TextFile tfc = new TextFile(superclasses, TextFile.R);
		tfc.readLine();
		String[] elemsf = tfc.readLineElems(TextFile.tab);
		while (elemsf != null) {
			superpop.put(elemsf[0], elemsf[2]);
			elemsf = tfc.readLineElems(TextFile.tab);
		}
		tfc.close();
		
		String nullpop = "NotDefined";
		int startk = 5;
		
		HashMap<String, HashSet<String>> samplesPerPop = new HashMap<>();
		HashMap<String, String> sampleToPop = new HashMap<>();
		TextFile tf = new TextFile(knownclasses, TextFile.R);
		
		tf.readLine();
		
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String sample = elems[0];
			String type = elems[1];
			type = superpop.get(type);
			HashSet<String> set = samplesPerPop.get(type);
			if (set == null) {
				set = new HashSet<>();
			}
			set.add(sample);
			samplesPerPop.put(type, set);
			sampleToPop.put(sample, type);
			
			elems = tf.readLineElems(TextFile.tab);
			
		}
		tf.close();
		
		System.out.println(sampleToPop.size() + " with known population.");
		System.out.println(samplesPerPop.size() + " populations total.");
		
		
		ArrayList<String> populations = new ArrayList<>();
		ArrayList<String> populationstmp = new ArrayList<>();
		populationstmp.addAll(samplesPerPop.keySet());
		Collections.sort(populationstmp);
		populations.add(nullpop);
		populations.addAll(populationstmp);
		
		
		HashMap<String, Integer> populationIndex = new HashMap<>();
		ArrayList<ArrayList<Double>> xvals = new ArrayList<>();
		ArrayList<ArrayList<Double>> yvals = new ArrayList<>();
		ArrayList<ArrayList<Double>> xvalsunlabeled = new ArrayList<>();
		ArrayList<ArrayList<Double>> yvalsunlabeled = new ArrayList<>();
		for (int i = 0; i < populations.size(); i++) {
			String pop = populations.get(i);
			populationIndex.put(pop, i);
			xvals.add(new ArrayList<>());
			xvalsunlabeled.add(new ArrayList<>());
			yvals.add(new ArrayList<>());
			yvalsunlabeled.add(new ArrayList<>());
		}
		
		ArrayList<Triple<String, Double, Double>> data = new ArrayList<Triple<String, Double, Double>>();
		ArrayList<Triple<String, Double, Double>> data2 = new ArrayList<Triple<String, Double, Double>>();
		TextFile tf2 = new TextFile(pcafile, TextFile.R);
		String ln = tf2.readLine();
		while (ln != null) {
			String pcadata = ln;
			while (pcadata.contains("  ")) {
				ln.replaceAll("  ", " ");
			}
			String[] pcaelems = pcadata.split(" ");
			String sample = pcaelems[0];
			Double pca1 = Double.parseDouble(pcaelems[2]);
			Double pca2 = Double.parseDouble(pcaelems[3]);
			Double pca3 = Double.parseDouble(pcaelems[4]);
			Double pca4 = Double.parseDouble(pcaelems[5]);
			
			data.add(new Triple<>(sample, pca1, pca2));
			data2.add(new Triple<>(sample, pca3, pca4));
			
			String population = sampleToPop.get(sample);
			if (population == null) {
				population = nullpop;
			}
			sampleToPop.put(sample, population);
			Integer popindex = populationIndex.get(population);
			if (population.equals(nullpop)) {
				xvalsunlabeled.get(popindex).add(pca1);
				yvalsunlabeled.get(popindex).add(pca2);
			} else {
				xvals.get(popindex).add(pca1);
				yvals.get(popindex).add(pca2);
			}
			ln = tf2.readLine();
		}
		tf2.close();
		
		DefaultTheme def = new DefaultTheme();
		Color[] colors = new Color[]{
				new Color(0, 0, 0),
				new Color(145, 145, 145),
				new Color(39, 121, 193),
				new Color(50, 189, 38),
				new Color(255, 211, 51),
				new Color(255, 51, 51),
		};
		def.setColors(colors);
		
		// scatterplot
		
		double[][] x = new double[populations.size()][];
		double[][] y = new double[populations.size()][];
		double[][] xunlab = new double[populations.size()][];
		double[][] yunlab = new double[populations.size()][];
		for (int i = 0; i < xvals.size(); i++) {
			x[i] = Primitives.toPrimitiveArr(xvals.get(i));
			y[i] = Primitives.toPrimitiveArr(yvals.get(i));
		}
		for (int i = 0; i < xvalsunlabeled.size(); i++) {
			xunlab[i] = Primitives.toPrimitiveArr(xvalsunlabeled.get(i));
			yunlab[i] = Primitives.toPrimitiveArr(yvalsunlabeled.get(i));
		}
		
		ScatterplotPanel p = new ScatterplotPanel(1, 1);
		p.setTitle("1000 genomes");
		p.setData(x, y);
		p.setDatasetLabels(populations.toArray(new String[0]));
		p.setLabels("PC1", "PC2");
		p.setPlotElems(true, true);
		p.setTheme(def);
		p.setAlpha(0.8f);
		ScatterplotPanel punlab = new ScatterplotPanel(1, 1);
		punlab.setTitle("Unlabeled samples");
		punlab.setData(xunlab, yunlab);
		punlab.setDatasetLabels(populations.toArray(new String[0]));
		punlab.setLabels("PC1", "PC2");
		punlab.setPlotElems(true, true);
		punlab.setTheme(def);
		punlab.setAlpha(0.8f);
		Grid grid = new Grid(500, 500, 1, 4, 100, 100);
		grid.addPanel(p);
		grid.addPanel(punlab);
		
		// knn population assignment
		HashMap<String, String> sampleToPopTmp = assignPopulation(data, sampleToPop, nullpop, startk);
		
		xvalsunlabeled = new ArrayList<>();
		yvalsunlabeled = new ArrayList<>();
		ArrayList<ArrayList<Double>> xvalsunlabeledpc3 = new ArrayList<>();
		ArrayList<ArrayList<Double>> yvalsunlabeledpc4 = new ArrayList<>();
		for (int i = 0; i < populations.size(); i++) {
			xvalsunlabeled.add(new ArrayList<>());
			yvalsunlabeled.add(new ArrayList<>());
			xvalsunlabeledpc3.add(new ArrayList<>());
			yvalsunlabeledpc4.add(new ArrayList<>());
		}
		
		TextFile assignmentout = new TextFile(output + "SampleAssignment.txt", TextFile.W);
		
		for (Triple<String, Double, Double> sample : data) {
			String prevpop = sampleToPop.get(sample.getLeft());
			if (prevpop.equals(nullpop)) {
				String pop = sampleToPopTmp.get(sample.getLeft());
				Integer popindex = populationIndex.get(pop);
				if (popindex == null) {
					System.out.println("Could not find population : " + pop + " for sample: " + sample.getLeft());
				} else {
					xvalsunlabeled.get(popindex).add(sample.getMiddle());
					yvalsunlabeled.get(popindex).add(sample.getRight());
					assignmentout.writeln(sample.getLeft() + "\t" + pop);
				}
			}
		}
		for (Triple<String, Double, Double> sample : data2) {
			String prevpop = sampleToPop.get(sample.getLeft());
			if (prevpop.equals(nullpop)) {
				String pop = sampleToPopTmp.get(sample.getLeft());
				Integer popindex = populationIndex.get(pop);
				if (popindex == null) {
					System.out.println("Could not find population : " + pop + " for sample: " + sample.getLeft());
				} else {
					xvalsunlabeledpc3.get(popindex).add(sample.getMiddle());
					yvalsunlabeledpc4.get(popindex).add(sample.getRight());
					assignmentout.writeln(sample.getLeft() + "\t" + pop);
				}
			}
		}
		assignmentout.close();
		xunlab = new double[populations.size()][];
		yunlab = new double[populations.size()][];
		double[][] xunlabpc3 = new double[populations.size()][];
		double[][] yunlabpc4 = new double[populations.size()][];
		for (int i = 0; i < xvalsunlabeled.size(); i++) {
			xunlab[i] = Primitives.toPrimitiveArr(xvalsunlabeled.get(i));
			yunlab[i] = Primitives.toPrimitiveArr(yvalsunlabeled.get(i));
			xunlabpc3[i] = Primitives.toPrimitiveArr(xvalsunlabeledpc3.get(i));
			yunlabpc4[i] = Primitives.toPrimitiveArr(yvalsunlabeledpc4.get(i));
		}
		
		Range datarange = punlab.getDataRange();
		
		punlab = new ScatterplotPanel(1, 1);
		punlab.setDataRange(datarange);
		punlab.setData(xunlab, yunlab);
		punlab.setTitle("Samples after labeling");
		punlab.setDatasetLabels(populations.toArray(new String[0]));
		punlab.setLabels("PC1", "PC2");
		punlab.setPlotElems(true, true);
		punlab.setTheme(def);
		punlab.setAlpha(0.8f);
		grid.addPanel(punlab);
		
		punlab = new ScatterplotPanel(1, 1);
		punlab.setDataRange(datarange);
		punlab.setData(xunlabpc3, yunlabpc4);
		punlab.setTitle("Samples after labeling");
		punlab.setDatasetLabels(populations.toArray(new String[0]));
		punlab.setLabels("PC3", "PC4");
		punlab.setPlotElems(true, true);
		punlab.setTheme(def);
		punlab.setAlpha(0.8f);
		grid.addPanel(punlab);
		
		grid.draw(output + "SampleAssignment.pdf");
	}
	
	private HashMap<String, String> assignPopulation(ArrayList<Triple<String, Double, Double>> data, HashMap<String, String> sampleToPop, String nullpop, int startk) {
		int samplesWithoutPop = 0;
		HashMap<String, String> sampleToPopTmp = new HashMap<>();
		
		for (Triple<String, Double, Double> sample : data) {
			String pop = sampleToPop.get(sample.getLeft());
			
			
			if (pop.equals(nullpop)) {
				int k = startk;
				boolean resolved = false;
				while (!resolved) {
					String[] neighbors = new String[k];
					double[] distances = new double[k];
					double maxdist = -Double.MAX_VALUE;
					int n = 0;
					
					
					// this sample needs assignment
					// get startk nearest neighbors that do have an assignment
					for (Triple<String, Double, Double> sample2 : data) {
						String pop2 = sampleToPop.get(sample2.getLeft());
						if (!pop2.equals(nullpop)) {
							double distance = Math.sqrt(sq(sample.getMiddle() - sample2.getMiddle()) + sq(sample.getRight() - sample2.getRight()));
							
							if (n < k) {
								neighbors[n] = sample2.getLeft();
								distances[n] = distance;
								if (distance > maxdist) {
									maxdist = distance;
								}
								n++;
							} else if (distance < maxdist) {
								boolean update = false;
								double newmax = -Double.MAX_VALUE;
								
								for (int i = 0; i < neighbors.length; i++) {
									if (!update && (distances[i] > distance)) {
										distances[i] = distance;
										neighbors[i] = sample2.getLeft();
										update = true;
									}
									
									if (distances[i] > newmax) {
										newmax = distances[i];
									}
								}
								maxdist = newmax;
							}
						}
					}
					
					HashMap<String, Integer> ctr = new HashMap<>();
					HashMap<String, Double> sumd = new HashMap<>();
					
					for (int i = 0; i < neighbors.length; i++) {
						String pop2 = sampleToPop.get(neighbors[i]);
						Integer ct = ctr.get(pop2);
						Double sum = sumd.get(pop2);
						if (ct == null) {
							ct = 1;
							sum = distances[i];
						} else {
							ct++;
							sum += distances[i];
						}
						ctr.put(pop2, ct);
						sumd.put(pop2, sum);
					}
					
					
					
					int maxct = 0;
					String maxpop = null;
					
					for (String key : ctr.keySet()) {
						Integer ct = ctr.get(key);
						if (ct > maxct) {
							maxpop = key;
							maxct = ct;
						} else if (ct == maxct) {
							// break ties using absolute sum of distance
							double sum = sumd.get(key);
							if (sum < sumd.get(maxpop)) {
								maxpop = key;
								maxct = ct;
							}
						}
					}
					
					// assign
					sampleToPopTmp.put(sample.getLeft(), maxpop);
					resolved = true;
				}
			} else {
				
				sampleToPopTmp.put(sample.getLeft(), pop);
			}
		}
		
		return sampleToPopTmp;
	}
	
	private double sq(double v) {
		return v * v;
	}
	
	
}
