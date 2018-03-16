/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package dashapca;

import java.io.IOException;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.ExpressionDataset;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.PCA;

/**
 *
 * @author harmjan
 */
public class DashaPCA {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        if (args.length < 1) {
            System.out.println("Usage: DashaPCA.jar expressionfile");
            System.exit(-1);
        }
        try {
            DashaPCA d = new DashaPCA();
            d.run(args[0]);
        } catch (IOException e) {
            e.printStackTrace();;
        }
        System.exit(-1);
    }

    public void run(String expressionFile) throws IOException {
        ExpressionDataset dataset = new ExpressionDataset(expressionFile);
        double[][] rawData = dataset.getRawData();

        System.out.print("- Calculating correlations between all " + dataset.getNrSamples() + " samples: ");
        double[][] correlationMatrix = new double[dataset.getNrSamples()][dataset.getNrSamples()];
        double probeCountMinusOne = dataset.getNrProbes() - 1;

        ProgressBar pb = new ProgressBar(dataset.getNrSamples());
        for (int f = 0; f < dataset.getNrSamples(); f++) {
            for (int g = f; g < dataset.getNrSamples(); g++) {
                double covarianceInterim = 0;
                for (int p = 0; p < dataset.getNrProbes(); p++) {
                    covarianceInterim += dataset.getRawData()[p][f] * dataset.getRawData()[p][g];
                }
                double covariance = covarianceInterim / probeCountMinusOne;
                correlationMatrix[f][g] = covariance;
                correlationMatrix[g][f] = covariance;
            }
            pb.iterate();
        }
        pb.close();

        System.out.println("- Performing PCA over samples:");
        Jama.EigenvalueDecomposition eig = PCA.eigenValueDecomposition(correlationMatrix);
        ExpressionDataset datasetEV = new ExpressionDataset(dataset.getNrSamples(), dataset.getNrSamples());
        datasetEV.setProbeNames(dataset.getSampleNames());
        double[] eigenValues = eig.getRealEigenvalues();
        System.out.println("Eigenvalue results:");

        TextFile out = new TextFile(expressionFile + ".PCAOverSamplesEigenvalues.txt.gz", TextFile.W);


        double cumExpVarPCA = 0;
        for (int pca = 0; pca < dataset.getNrSamples(); pca++) {
            double expVarPCA = PCA.getEigenValueVar(eigenValues, pca);
            double[] pca1ExpEigenVector = PCA.getEigenVector(eig, eigenValues, pca);
            for (int s = 0; s < dataset.getNrSamples(); s++) {
                datasetEV.getRawData()[s][pca] = pca1ExpEigenVector[s];
            }
            int pcaNr = pca + 1;
            cumExpVarPCA += expVarPCA;
            out.write(pcaNr + "\t" + expVarPCA + "\t" + cumExpVarPCA + "\n");
            datasetEV.getSampleNames()[pca] = "Comp" + String.valueOf(pcaNr);
            System.out.println("PCA over samples:\t" + pcaNr + "\t" + eigenValues[eigenValues.length - 1 - pca] + "\t" + expVarPCA);
        }
        out.close();

        datasetEV.save(expressionFile + ".PCAOverSamplesEigenvectors.txt.gz", true);

        ExpressionDataset datasetEVT = new ExpressionDataset(dataset.getNrSamples(), dataset.getNrSamples());
        datasetEVT.setProbeNames(datasetEV.getProbeNames());
        datasetEVT.setSampleNames(datasetEV.getSampleNames());
        datasetEVT.setRawData(datasetEV.getRawData());
        datasetEVT.transposeDataset();

        datasetEVT.save(expressionFile + ".PCAOverSamplesEigenvectorsTransposed.txt.gz", true);
        datasetEVT = null;


        System.out.println("Calculating PCs:");
        System.out.println(" - Initializing PCA matrix:");
        ExpressionDataset datasetPCAOverSamplesPCAs = new ExpressionDataset(dataset.getNrProbes(), dataset.getNrSamples());
        datasetPCAOverSamplesPCAs.setProbeNames(dataset.getProbeNames());
        for (int s = 0; s < dataset.getNrSamples(); s++) {
            datasetPCAOverSamplesPCAs.getSampleNames()[s] = "Comp" + String.valueOf(s + 1);
        }
        for (int p = 0; p < dataset.getNrProbes(); p++) {
            for (int t = 0; t < dataset.getNrSamples(); t++) {
                datasetPCAOverSamplesPCAs.getRawData()[p][t] = 0;
            }
        }

        System.out.println(" - Calculating per probe the PCA scores:");
        pb = new ProgressBar(dataset.getNrProbes());
        for (int probe = 0; probe < dataset.getNrProbes(); probe++) {
            for (int sample1 = 0; sample1 < dataset.getNrSamples(); sample1++) {
                for (int sample2 = 0; sample2 < dataset.getNrSamples(); sample2++) {
                    double probeCoefficient = datasetEV.getRawData()[sample2][sample1];
                    datasetPCAOverSamplesPCAs.getRawData()[probe][sample1] += probeCoefficient * dataset.getRawData()[probe][sample2];
                }
            }
            pb.iterate();
        }
        pb.close();

        System.out.println(" - Saving PCA scores:");
        datasetPCAOverSamplesPCAs.save(expressionFile + ".PCAOverSamplesPrincipalComponents.txt.gz", true);

        System.out.println("Done\n");
    }
}
