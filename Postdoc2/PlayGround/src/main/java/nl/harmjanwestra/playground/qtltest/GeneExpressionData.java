package nl.harmjanwestra.playground.qtltest;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class GeneExpressionData {

	public double[][] data;
	public String[] genes;
	public String[] samples;
	public HashMap<String, Integer> sampleMap;
	public HashMap<String, Integer> geneMap;

	public GeneExpressionData(String geneExpressionDataFile, Set<String> geneSelection, Set<String> requestedSamples) throws IOException {
		System.out.println("Loading expression data from: " + geneExpressionDataFile);
		if (requestedSamples != null) {
			System.out.println("Max number of samples: " + requestedSamples.size());
		}
		if (geneSelection != null) {
			System.out.println("Max number of genes: " + geneSelection.size());
		}

		TextFile tf = new TextFile(geneExpressionDataFile, TextFile.R);
		String[] header = tf.readLineElems(TextFile.tab);
		boolean[] includeColumn = new boolean[header.length];

		ArrayList<String> sampleTmp = new ArrayList<>();
		for (int i = 1; i < header.length; i++) {
			String sample = header[i];
			if (requestedSamples == null || requestedSamples.contains(sample)) {
				includeColumn[i] = true;
				sampleTmp.add(sample);
			}
		}
		System.out.println(sampleTmp.size() + " samples found.");
		samples = sampleTmp.toArray(new String[0]);
		sampleMap = Util.hash(sampleTmp);

		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<double[]> dataList = new ArrayList<>();
		ArrayList<String> genetmp = new ArrayList<>();
		int lctr = 0;
		while (elems != null) {
			String gene = Strings.cache(elems[0]);
			if (geneSelection == null || geneSelection.contains(gene)) {
				double[] dataln = new double[samples.length];
				genetmp.add(gene);
				int sctr = 0;
				for (int i = 1; i < elems.length; i++) {
					if (includeColumn[i]) {
						double d = Double.parseDouble(elems[i]);
						dataln[sctr] = d;
						sctr++;
					}
				}
				dataList.add(dataln);
			}
			lctr++;
			if (lctr % 2000 == 0) {
				System.out.print(lctr + " lines parsed, " + dataList.size() + " genes loaded.\r");
//                break;
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(lctr + " lines parsed, " + dataList.size() + " genes loaded.");

		genes = genetmp.toArray(new String[0]);
		data = new double[genes.length][0];
		for (int g = 0; g < genes.length; g++) {
			data[g] = dataList.get(g);
		}
		geneMap = Util.hash(genetmp);
	}

}
