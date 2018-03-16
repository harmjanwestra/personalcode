/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

/**
 *
 * @author harmjan
 */
public class TriTyperGenotypeCompare {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        TriTyperGenotypeCompare fd = new TriTyperGenotypeCompare();
        String ds1 = "/Volumes/iSnackHD/TonuData/1000G_Nov2011_release/TriTyper/";
        String ds2 = "/Volumes/iSnackHD/TonuData/1000G_Nov2011_release/TriTyper3/";
        try {
//            fd.test(ds1);
//            fd.test(ds2);
            fd.run(ds1, ds2);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void test(String ds1loc) throws IOException {
        TriTyperGenotypeData ds1 = new TriTyperGenotypeData(ds1loc);
        SNPLoader ds1loader = ds1.createSNPLoader();
        String[] snps1 = ds1.getSNPs();
        
        int nullctr = 0;
        for(int snp1Id=0; snp1Id<snps1.length; snp1Id++){
            SNP snpObj1 = ds1.getSNPObject(snp1Id);
            ds1loader.loadGenotypes(snpObj1);
            double cr = snpObj1.getCR();
            if(cr == 0d){
                nullctr++;
            }
            
        }
        
        System.out.println(nullctr+" out of "+snps1.length + " have zero callrate?");
        
        
    }
    
    public void run(String ds1loc, String ds2loc) throws IOException {
        TriTyperGenotypeData ds1 = new TriTyperGenotypeData(ds1loc);
        TriTyperGenotypeData ds2 = new TriTyperGenotypeData(ds2loc);

        SNPLoader ds1loader = ds1.createSNPLoader();
        SNPLoader ds2loader = ds2.createSNPLoader();

        String[] snps1 = ds1.getSNPs();
        String[] snps2 = ds2.getSNPs();

        String[] individuals1 = ds1.getIndividuals();
        int sharedInds = 0;
        int[] toInd2 = new int[individuals1.length];
        for (int ind1 = 0; ind1 < individuals1.length; ind1++) {
            Integer ind2Id = ds2.getIndividualToId().get(individuals1[ind1]);
            if (ind2Id != null) {
                toInd2[ind1] = ind2Id;
                sharedInds++;
            } else {
                toInd2[ind1] = -1;
            }
        }

        System.out.println(sharedInds + " shared inds");



        double correlationSum = 0;
        int sharedSNPs = 0;
        ProgressBar pb = new ProgressBar(snps1.length);
        for (int snp1Id = 0; snp1Id < snps1.length; snp1Id++) {
            String snp = snps1[snp1Id];
            Integer snp2Id = ds2.getSnpToSNPId().get(snp);
            if (snp2Id != null) {
                sharedSNPs++;
                SNP snpObj1 = ds1.getSNPObject(snp1Id);
                SNP snpObj2 = ds2.getSNPObject(snp2Id);

                ds1loader.loadGenotypes(snpObj1);
                ds2loader.loadGenotypes(snpObj2);

                byte[] geno1 = snpObj1.getGenotypes();
                byte[] geno2 = snpObj2.getGenotypes();

                int ctr = 0;
                double[] geno1dbl = new double[sharedInds];
                double[] geno2dbl = new double[sharedInds];
                for (int ind1Id = 0; ind1Id < individuals1.length; ind1Id++) {
                    int ind2Id = toInd2[ind1Id];
                    if (ind2Id != -1) {
                        geno1dbl[ctr] = geno1[ind1Id];
                        geno2dbl[ctr] = geno2[ind2Id];
                        ctr++;
                    }
                }
                snpObj1.clearGenotypes();
                snpObj2.clearGenotypes();

                double corr = Math.abs(JSci.maths.ArrayMath.correlation(geno1dbl, geno2dbl));
                if (corr == 0.0) {
                    ctr = 0;
                    for (int ind1Id = 0; ind1Id < individuals1.length; ind1Id++) {
                        int ind2Id = toInd2[ind1Id];
                        if (ind2Id != -1) {
                            System.out.println(individuals1[ind1Id]+"\t"+geno1dbl[ctr]+"\t"+geno2dbl[ctr]);
                            
                            ctr++;
                        }
                    }
                    System.out.println();
                    System.exit(0);
                }
                correlationSum += corr;
            }
            pb.iterate();
        }
        pb.close();



        System.out.println("Shared SNPs: " + sharedSNPs);
        System.out.println("Average absolute correlation: " + (correlationSum / sharedSNPs));



    }
}
