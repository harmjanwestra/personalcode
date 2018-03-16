/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package causalitytest;

import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Regression;

/**
 *
 * @author harmjan
 */
public class CausalityTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        CausalityTest.test();
    }

    public static void test() {
        // model genotypes with 50% MAF
        int nrgenotypes = 1240;
        double[] genotypes = new double[nrgenotypes];

        for (int s=0; s<nrgenotypes; s++) {
            genotypes[s] = 0;
            if (Math.random()>0.5) {
                genotypes[s]++;
            }
            if (Math.random()>0.5) {
                genotypes[s]++;
            }
        }
        

        // now model the cis effect
        double[] cisProbe = new double[nrgenotypes];
        for (int i = 0; i < nrgenotypes; i++) {
            cisProbe[i] = genotypes[i] + 25d * (Math.random() * Math.random() * Math.random() * Math.random() * Math.random());
        }

        // now model the independent trans effect
        double[] independentTransProbe = new double[nrgenotypes];
        for (int i = 0; i < nrgenotypes; i++) {
            independentTransProbe[i] = genotypes[i] + 25d * (Math.random() * Math.random() * Math.random() * Math.random() * Math.random());
        }


        // now model the dependent trans effect
        double[] dependentTransProbe = new double[nrgenotypes];
        for (int i = 0; i < nrgenotypes; i++) {
            dependentTransProbe[i] = cisProbe[i] + 25d * (Math.random() * Math.random() * Math.random() * Math.random() * Math.random());
        }

        // now calculate the different statistics

        // regress cis against dependent effect and independent trans
        double[] snpCisEffect = Regression.getLinearRegressionCoefficients(genotypes, cisProbe);
        double[] snpDependentTransEffect = Regression.getLinearRegressionCoefficients(genotypes, dependentTransProbe);
        double[] snpInDependentTransEffect = Regression.getLinearRegressionCoefficients(genotypes, independentTransProbe);

        double[] cisDependentTransEffect = Regression.getLinearRegressionCoefficients(dependentTransProbe, cisProbe);
        double[] cisInDependentTransEffect = Regression.getLinearRegressionCoefficients(independentTransProbe, cisProbe);

        System.out.println("Cis effect beta:\t\t" + snpCisEffect[0]+"\t"+Correlation.correlate(genotypes, cisProbe));
        System.out.println("Dependent trans effect beta:\t" + snpDependentTransEffect[0]+"\t"+Correlation.correlate(genotypes, dependentTransProbe));
        System.out.println("InDependent trans effect beta:\t" + snpInDependentTransEffect[0]+"\t"+Correlation.correlate(genotypes, independentTransProbe));
        

        // now remove the snp effect from both the cis and the trans effects
        double[] resCis = new double[nrgenotypes];
        double[] resTransDependent = new double[nrgenotypes];
        double[] resTransIndependent = new double[nrgenotypes];

        for (int i = 0; i < nrgenotypes; i++) {
            resCis[i] = cisProbe[i] - snpCisEffect[0] * genotypes[i];
            resTransDependent[i] = dependentTransProbe[i] - snpDependentTransEffect[0] * genotypes[i];
            resTransIndependent[i] = independentTransProbe[i] - snpInDependentTransEffect[0] * genotypes[i];
        }

        System.out.println("\n\n");
        
        double corrCisTransDependent = JSci.maths.ArrayMath.correlation(cisProbe, dependentTransProbe);
        double corrCisTransIndependent = JSci.maths.ArrayMath.correlation(cisProbe, independentTransProbe);
        double corrResCisResTransDependent = JSci.maths.ArrayMath.correlation(resCis, resTransDependent);
        double corrResCisResTransIndependent = JSci.maths.ArrayMath.correlation(resCis, resTransIndependent);
        System.out.println("corrCisTransDependent:\t"+ corrCisTransDependent);
        System.out.println("corrResCisResTransDependent:\t"+ corrResCisResTransDependent);

        System.out.println("\n\n");
        
        System.out.println("corrCisTransIndependent:\t"+ corrCisTransIndependent);
        System.out.println("corrResCisResTransIndependent:\t"+ corrResCisResTransIndependent);
        
        
        double cisInDependentTransCorrelation = Correlation.correlate(cisProbe, independentTransProbe);
        double cisDependentTransCorrelation = Correlation.correlate(cisProbe, dependentTransProbe);

        double cisInDependentTransCorrelationAfterRegression = Correlation.correlate(resCis, resTransIndependent);
        double cisDependentTransCorrelationAfterRegression = Correlation.correlate(resCis, resTransDependent);


        
        System.out.println("\n\n");
        
        System.out.println("Before regression: \t\tdep: " + cisDependentTransCorrelation + "\tindep: " + cisInDependentTransCorrelation);
        System.out.println("After regression: \t\tdep: " + cisDependentTransCorrelationAfterRegression + "\tindep: " + cisInDependentTransCorrelationAfterRegression);

    }

    public static double getRandomGenotype() {
        double genotype = Math.random() * 2;
        double iGenotype = 0;
        if (genotype < 0.5) {
            iGenotype = 0;
        } else if (genotype > 1.5) {
            iGenotype = 2;
        } else {
            iGenotype = 1;
        }
        return iGenotype;
    }
}
