package nl.harmjanwestra.playground.biogen.freeze2dot1.interactions;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

public class InteractionSim {


    public static void main(String[] args) {
        InteractionSim s = new InteractionSim();
        s.run(1000, 0.10, -0.5, 0.7, 0.3);
    }

    public void run(int nrGenotypes, double maf, double snpfx, double covfx, double intfx) {

        double[] genotypes = new double[nrGenotypes];

        double nrA = maf * nrGenotypes;
        double nrhets = 0.75 * nrA;

        double nrB = (1 - maf) * nrGenotypes;

        int set = 0;
        for (int i = 0; i < nrA - nrhets; i++) {
            genotypes[i] = 0;
            set++;
        }
        for (int i = set; i < (nrA - nrhets) + nrhets; i++) {
            genotypes[i] = 1;
            set++;
        }
        for (int i = set; i < nrGenotypes; i++) {
            genotypes[i] = 2;
            set++;
        }

        int ctrA = 0;
        int ctrB = 0;
        for (int i = 0; i < nrGenotypes; i++) {
            if (genotypes[i] == 0) {
                ctrA += 2;
            } else if (genotypes[i] == 1) {
                ctrA += 1;
                ctrB += 1;
            } else {
                ctrB += 2;
            }
        }

        double[] cov = new double[nrGenotypes];
        for (int i = 0; i < cov.length; i++) {
            cov[i] = Math.random();
        }

        System.out.println(nrGenotypes + "\t" + ctrA + "\t" + ctrB + "\t" + ((double) ctrA / (nrGenotypes * 2)));
        double[] expression = new double[nrGenotypes];
        double[] noise = new double[nrGenotypes];

        for (int i = 0; i < expression.length; i++) {
            noise[i] = Math.random() + Math.random() + Math.random();
            expression[i] = (snpfx * genotypes[i]) + (covfx * cov[i]) + (intfx * genotypes[i] * cov[i]) + noise[i];
        }

        double[][] x = new double[nrGenotypes][3];
        System.out.println("Input:");
        System.out.println("SNPFX: " + snpfx);
        System.out.println("CovFX: " + covfx);
        System.out.println("IntFx: " + intfx);
        for (int i = 0; i < x.length; i++) {
            x[i][0] = genotypes[i];
            x[i][1] = cov[i];
            x[i][2] = cov[i] * genotypes[i];
        }

        OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
        ols.newSampleData(expression, x);
        double[] betas = ols.estimateRegressionParameters();
        for (int d = 0; d < betas.length; d++) {
            if (d == 0) {
                System.out.println(d + "\tIntercept: " + betas[d]);
            } else if (d == 1) {
                System.out.println(d + "\tSNP: " + betas[d]);
            } else if (d == 2) {
                System.out.println(d + "\tCov: " + betas[d]);
            } else {
                System.out.println(d + "\tCov*SNP: " + betas[d]);
            }

        }
        System.out.println("Rsquared: " + ols.calculateAdjustedRSquared());


        for (int i = 0; i < expression.length; i++) {
            genotypes[i] = 2 - genotypes[i];
//            expression[i] = (snpfx * genotypes[i]) + (covfx *  cov[i]) + (intfx * genotypes[i] * cov[i]); // + Math.random() + Math.random() + Math.random();
        }

        for (int i = 0; i < x.length; i++) {
            x[i][0] = genotypes[i];
            x[i][1] = cov[i];
            x[i][2] = cov[i] * genotypes[i];
        }
        System.out.println();
        System.out.println("Alleles flipped:");
        ols = new OLSMultipleLinearRegression();
        ols.newSampleData(expression, x);
        betas = ols.estimateRegressionParameters();
        for (int d = 0; d < betas.length; d++) {
            if (d == 0) {
                System.out.println(d + "\tIntercept: " + betas[d]);
            } else if (d == 1) {
                System.out.println(d + "\tSNP: " + betas[d]);
            } else if (d == 2) {
                System.out.println(d + "\tCov: " + betas[d]);
            } else {
                System.out.println(d + "\tCov*SNP: " + betas[d]);
            }
        }
        System.out.println("Rsquared: " + ols.calculateAdjustedRSquared());

        for (int i = 0; i < expression.length; i++) {
            expression[i] = (snpfx * genotypes[i]) + (covfx * cov[i]) + (intfx * genotypes[i] * cov[i]) + noise[i];
        }

        for (int i = 0; i < x.length; i++) {
            x[i][0] = genotypes[i];
            x[i][1] = cov[i];
            x[i][2] = cov[i] * genotypes[i];
        }

        System.out.println();
        System.out.println("Alleles flipped, expression recalculated:");
        ols = new OLSMultipleLinearRegression();
        ols.newSampleData(expression, x);
        betas = ols.estimateRegressionParameters();
        for (int d = 0; d < betas.length; d++) {
            if (d == 0) {
                System.out.println(d + "\tIntercept: " + betas[d]);
            } else if (d == 1) {
                System.out.println(d + "\tSNP: " + betas[d]);
            } else if (d == 2) {
                System.out.println(d + "\tCov: " + betas[d]);
            } else {
                System.out.println(d + "\tCov*SNP: " + betas[d]);
            }
        }
        System.out.println("Rsquared: " + ols.calculateAdjustedRSquared());

    }
}
