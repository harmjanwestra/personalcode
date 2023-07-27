package nl.harmjanwestra.playground.coeqtl;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

public class CoEQTLSizeTable {


	public static void main(String[] args) {
		int sampleStepSize = 50;
		int maxSampleSize = 1000 + sampleStepSize;
		int geneStepSize = 500;
		int maxGenes = 20000;
		int precision = 8; // 8 byte float
		String header = "NrGenes\tNrPairs";
		for (int sampleSize = sampleStepSize; sampleSize < maxSampleSize; sampleSize += sampleStepSize) {
			header += "\tn" + sampleSize;
		}
//		System.out.println(header);
//
//
//		for (int nrGenes = geneStepSize; nrGenes < maxGenes; nrGenes += geneStepSize) {
//			long pairs = ((nrGenes * nrGenes) / 2) - nrGenes;
//			String outln = nrGenes + "\t" + pairs;
//			for (int sampleSize = sampleStepSize; sampleSize < maxSampleSize; sampleSize += sampleStepSize) {
//				long values = pairs * sampleSize;
//				long bytes = values * precision;
//				outln += "\t" + Gpio.humanizeFileSize(bytes);
//			}
//			System.out.println(outln);
//		}

		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println();
		try {
			TextFile tf = new TextFile("D:\\Sync\\SyncThing\\Postdoc2\\2023-coEQTL\\Datasets.txt", TextFile.R);

			ArrayList<Pair<String, Integer>> datasets = new ArrayList<Pair<String, Integer>>();
			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {

				Pair<String, Integer> p = new Pair<>(elems[0], Integer.parseInt(elems[1]));
				datasets.add(p);

				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			header = "NrGenes\tNrPairs";

			int sumN = 0;
			for (Pair<String, Integer> dataset : datasets) {
				header += "\t" + dataset.getLeft() + " (n=" + dataset.getRight() + ")";
				sumN += dataset.getRight();
			}
//			header += "\tMeta-Analysis (n=" + sumNrela + ")";

			System.out.println(header);
			for (int nrGenes = geneStepSize; nrGenes < maxGenes; nrGenes += geneStepSize) {
				long pairs = ((nrGenes * nrGenes) / 2) - nrGenes;
				String outln = nrGenes + "\t" + pairs;
				long sum = 0;
				for (Pair<String, Integer> dataset : datasets) {
					int sampleSize = dataset.getRight();
					long values = pairs * sampleSize;
					long bytes = values * precision;
					sum += bytes;
					outln += "\t" + Gpio.humanizeFileSize(bytes);
				}
				outln += "\t" + Gpio.humanizeFileSize(sum);
				System.out.println(outln);
			}

		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
