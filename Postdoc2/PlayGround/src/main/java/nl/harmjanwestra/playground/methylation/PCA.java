package nl.harmjanwestra.playground.methylation;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.PCAojAlgo;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.ArrayList;

public class PCA {

	public void run(String in, int nrOfPCsToCalculate) throws Exception {

		DoubleMatrixDataset<String, String> dataset = DoubleMatrixDataset.loadDoubleData(in);
		PCAojAlgo pcaObj = new PCAojAlgo();

		pcaObj.eigenValueDecomposition(dataset.getMatrix().toArray());


		DoubleMatrixDataset<String, String> datasetEV = new DoubleMatrixDataset<String, String>(dataset.columns(), nrOfPCsToCalculate);
		datasetEV.setRowObjects(dataset.getColObjects());
		datasetEV.setColObjects(new ArrayList<>());
		double[] eigenValues = pcaObj.getRealEigenValues();
		System.out.println("Eigenvalue results:");

		System.out.println("PCA\tPCANr\tEigenValue\tExplainedVariance\tTotalExplainedVariance");

		TextFile out = new TextFile(in + ".PCAOverSamplesEigenvalues.txt.gz", TextFile.W);
		double cumExpVarPCA = 0;

		out.writeln("PCA\tPCANr\tEigenValue\tExplainedVariance\tTotalExplainedVariance");

		ArrayList<String> evcolnames = new ArrayList<>();
		for (int pca = 0; pca < nrOfPCsToCalculate; pca++) {
			double expVarPCA = pcaObj.getEigenValueVar(pca);

			double[] pca1ExpEigenVector = pcaObj.getEigenVector(pca);
			for (int s = 0; s < dataset.columns(); s++) {
				datasetEV.setElementQuick(s, pca, pca1ExpEigenVector[s]);
			}

			int pcaNr = pca + 1;
			cumExpVarPCA += expVarPCA;
			out.write(pcaNr + "\t" + eigenValues[pca] + "\t" + expVarPCA + "\t" + cumExpVarPCA + "\n");
			evcolnames.add(pca, "Comp" + pcaNr);
			if (pca < 10) {
				System.out.println("PCA:\t" + pcaNr + "\t" + eigenValues[pca] + "\t" + expVarPCA + "\t" + cumExpVarPCA);
			} else if (pca == 10) {
				System.out.println("Remaining eigenvalues in this file: " + in + ".PCAOverSamplesEigenvalues.txt.gz");
			}
		}

		datasetEV.setColObjects(evcolnames);

		datasetEV.save(in + ".PCAOverSamplesEigenvectors.txt.gz");


	}
}
