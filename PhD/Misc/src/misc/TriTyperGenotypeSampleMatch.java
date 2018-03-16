/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.Arrays;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.CompareAllelicDirections;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class TriTyperGenotypeSampleMatch {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        String ds1loc = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/Hap2ImputedGenotypes/";
        String ds2loc = "/Volumes/iSnackHD/TonuData/TT/";
        String outdir = "/Volumes/iSnackHD/TonuData/";
        
        try {
            TriTyperGenotypeData ds1 = new TriTyperGenotypeData(ds1loc);

            TriTyperGenotypeData ds2 = new TriTyperGenotypeData(ds2loc);

            SNPLoader loader1 = ds1.createSNPLoader();
            SNPLoader loader2 = ds2.createSNPLoader();



            int[][] matrix = new int[ds1.getIndividuals().length][ds2.getIndividuals().length];
            // find snps in both ds with a high MAF

            String[] individuals1 = ds1.getIndividuals();
            int nrInds1 = individuals1.length;
            String[] individuals2 = ds2.getIndividuals();
            int nrInds2 = individuals2.length;

            CompareAllelicDirections d = new CompareAllelicDirections();

            String[] snps1 = ds1.getSNPs();
            int nrSNPs = 0;
            for (int s = 0; s < 100000; s++) {
                Integer s2 = ds2.getSnpToSNPId().get(snps1[s]);
                if (s2 != null) {

                    SNP snpobj = ds1.getSNPObject(s);
                    loader1.loadGenotypes(snpobj);
                    if (snpobj.getMAF() > 0.25 && snpobj.getHWEP() > 0.00001) {
                        SNP snpobj2 = ds2.getSNPObject(s2);
                        loader2.loadGenotypes(snpobj2);

                        byte[] g1 = snpobj.getGenotypes();
                        byte[] g2 = snpobj2.getGenotypes();


                        Boolean[] flip = CompareAllelicDirections.compare(new SNP[]{snpobj, snpobj2});

                        if (flip!=null && flip[1] != null) {
                            nrSNPs++;

                            if (nrSNPs % 1000 == 0) {
                                System.out.println(nrSNPs + " snps parsed... out of " + s);
                            }
                            for (int i = 0; i < nrInds1; i++) {
                                for (int j = 0; j < nrInds2; j++) {
                                    if (!flip[1] && g1[i] == g2[j]) {
                                        // do true stuff, no flip genotypes equal
                                        matrix[i][j]++;
                                    } else if (g1[i] == 1 && g2[j] == 1) {
                                        // heterozygotes
                                        matrix[i][j]++;
                                    } else if (flip[1] && ((g1[i] == 0 && g2[j] == 2) || (g1[i] == 2 && g2[j] == 1))) {
                                        matrix[i][j]++;
                                    }
                                }
                            }
                        }
                        snpobj2.clearGenotypes();
                    }
                    snpobj.clearGenotypes();
                }

            }

            double[][] perc = new double[nrInds1][nrInds2];
            for(int i=0; i<nrInds1; i++){
                for(int j=0; j<nrInds2;j++){
                    perc[i][j] = (double) matrix[i][j] / nrSNPs;
                }
                
            }
            
            DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<String, String>(perc, Arrays.asList(individuals1), Arrays.asList(individuals2));
            dsout.save(outdir+"IBS.txt");
            for (int i = 0; i < nrInds1; i++) {

                double max = -1;
                int nr = -1;
                for (int j = 0; j < nrInds2; j++) {
                    if (perc[i][j] > max) {
                        nr = j;
                        max = perc[i][j];
                    }
                }
                System.out.println(individuals1[i] + "\t" + individuals2[nr] + "\t" + matrix[i][nr] + "\t" + ((double) matrix[i][nr] / nrSNPs));
            }
            
            


        } catch (IOException e) {

            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
