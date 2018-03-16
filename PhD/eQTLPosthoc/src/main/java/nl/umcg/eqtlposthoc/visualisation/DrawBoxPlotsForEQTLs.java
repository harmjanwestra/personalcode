/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.visualisation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.graphics.ViolinBoxPlot;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Descriptives;

/**
 *
 * @author harmjan
 */
public class DrawBoxPlotsForEQTLs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        TriTyperGeneticalGenomicsDatasetSettings[] settings = new TriTyperGeneticalGenomicsDatasetSettings[2];

        settings[0] = new TriTyperGeneticalGenomicsDatasetSettings();
        settings[0].cisAnalysis = true;
        settings[0].transAnalysis = true;
        settings[0].expressionLocation = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/ExpressionData/ExpressionData.txt.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt";
        settings[0].expressionplatform = "HT12v3";
        settings[0].genotypeLocation = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/Hap2ImputedGenotypes/";
        settings[0].genotypeToExpressionCoupling = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/Hap2ImputedGenotypes/GenotypeExpressionCoupling.txt";
        settings[0].logtransform = false;
        settings[0].name = "EGCUT";
        settings[0].probeannotation = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-HT12v3.txt";
        settings[0].quantilenormalize = false;
        settings[0].tsProbesConfine = null;
//
        settings[1] = new TriTyperGeneticalGenomicsDatasetSettings();
        settings[1].cisAnalysis = true;
        settings[1].transAnalysis = true;
        settings[1].expressionLocation = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData4SNPPCRemoved/ExpressionData.txt.gz.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt";
        settings[1].expressionplatform = "HT12v3";
        settings[1].genotypeLocation = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/";
        settings[1].genotypeToExpressionCoupling = "";
        settings[1].logtransform = false;
        settings[1].name = "Fehrmann HT12v3";
        settings[1].probeannotation = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-HT12v3.txt";
        settings[1].quantilenormalize = false;
        settings[1].tsProbesConfine = null;


//        settings[2] = new TriTyperGeneticalGenomicsDatasetSettings();
//        settings[2].cisAnalysis = true;
//        settings[2].transAnalysis = true;
//        settings[2].expressionLocation = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodH8v2/BloodH8v2OriginalExpressionDataCorrectedFor4GWASPCs/ExpressionData.txt.gz.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz";
//        settings[2].expressionplatform = "HT12v3";
//        settings[2].genotypeLocation = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodH8v2/BloodH8v2ImputeTriTyper/";
//        settings[2].genotypeToExpressionCoupling = "";
//        settings[2].logtransform = false;
//        settings[2].name = "Fehrmann H8v2";
//        settings[2].probeannotation = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-HT12v3.txt";
//        settings[2].quantilenormalize = false;
//        settings[2].tsProbesConfine = null;


        String[] snps = new String[]{"rs921320"};
        String[] probes = new String[]{"3360093"};


        DrawBoxPlotsForEQTLs d = new DrawBoxPlotsForEQTLs();
        try {
            d.run(settings, snps, probes, "/Volumes/iSnackHD/Data/Projects/SuzanneVanSommeren/eQTLResults/2013-04-23-COMMD1Locus-WoH8v2/plots/");
        } catch (IOException ex) {
            Logger.getLogger(DrawBoxPlotsForEQTLs.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(DrawBoxPlotsForEQTLs.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public double[][] centerAndScale(double[][] rawData) throws IOException {

        System.out.println("Standardizing probe mean and standard deviation");
        for (int p = 0; p < rawData.length; p++) {
            double mean = Descriptives.mean(rawData[p]);
            double stdev = Math.sqrt(Descriptives.variance(rawData[p], mean));
            for (int s = 0; s < rawData[p].length; s++) {
                rawData[p][s] -= mean;
                rawData[p][s] /= stdev;
            }
        }

        return rawData;
    }

    public void run(TriTyperGeneticalGenomicsDatasetSettings[] settings, String[] snps, String[] probes, String outdir) throws IOException, Exception {

        outdir = Gpio.formatAsDirectory(outdir);
        if (!Gpio.exists(outdir)) {
            Gpio.createDir(outdir);
        }

        eqtlmappingpipeline.normalization.Normalizer n = new eqtlmappingpipeline.normalization.Normalizer();
        // load datasets
        TriTyperGeneticalGenomicsDataset[] ds = new TriTyperGeneticalGenomicsDataset[settings.length];
        String[] datasetNames = new String[ds.length];
        for (int i = 0; i < ds.length; i++) {
            ds[i] = new TriTyperGeneticalGenomicsDataset(settings[i]);

            centerAndScale(ds[i].getExpressionData().getMatrix());

            datasetNames[i] = settings[i].name;
        }


        // 
        for (int s = 0; s < snps.length; s++) {

            String snp = snps[s];


            String allele1FirstSNP = null;
            String allele2FirstSNP = null;

            String minorAlleleFirstSNP = null;

            for (int p = 0; p < probes.length; p++) {
                String probe = probes[s];
                double[][][] vals = new double[ds.length][3][0];
                String[][] genotypeLabels = new String[ds.length][3];
                int dsctr = 0;
                for (int d = 0; d < ds.length; d++) {
                    // check whether the SNP is there
                    TriTyperGenotypeData gt = ds[d].getGenotypeData();
                    Integer snpId = gt.getSnpToSNPId().get(snp);
                    Integer probeId = ds[d].getExpressionData().getProbeToId().get(probe);
                    boolean flipAlleles = false;
                    boolean incompatible = false;
                    if (snpId == null || probeId == null) {
                        vals[d] = null;
                        System.out.println("SNP or Probe not in dataset: " + settings[d].name + "\t" + snp + "  - " + snpId + "\t" + probe + " - " + probeId);
                    } else {
                        dsctr++;
                        ArrayList<Double> valsgt0 = new ArrayList<Double>();
                        ArrayList<Double> valsgt1 = new ArrayList<Double>();
                        ArrayList<Double> valsgt2 = new ArrayList<Double>();


                        // collect values..

                        SNP snpObj = gt.getSNPObject(snpId);
                        SNPLoader loader = gt.createSNPLoader();
                        loader.loadGenotypes(snpObj);

                        if (allele1FirstSNP == null) {
                            flipAlleles = false;
                            allele1FirstSNP = BaseAnnot.toString(snpObj.getAlleles()[0]);
                            allele2FirstSNP = BaseAnnot.toString(snpObj.getAlleles()[1]);

                            minorAlleleFirstSNP = BaseAnnot.toString(snpObj.getMinorAllele());
                            genotypeLabels[d][0] = BaseAnnot.toString(snpObj.getAlleles()[0]) + BaseAnnot.toString(snpObj.getAlleles()[0]);
                            genotypeLabels[d][1] = BaseAnnot.toString(snpObj.getAlleles()[0]) + BaseAnnot.toString(snpObj.getAlleles()[1]);
                            genotypeLabels[d][2] = BaseAnnot.toString(snpObj.getAlleles()[1]) + BaseAnnot.toString(snpObj.getAlleles()[1]);

                        } else {

                            String allele1SNP = BaseAnnot.toString(snpObj.getAlleles()[0]);
                            String allele2SNP = BaseAnnot.toString(snpObj.getAlleles()[1]);
                            String minorAlleleSNP = BaseAnnot.toString(snpObj.getMinorAllele());

                            // what if A/T SNP?
                            if (allele1SNP.equals(BaseAnnot.getComplement(allele2SNP))) {
                                // compare minor alleles
                                if (minorAlleleSNP.equals(minorAlleleFirstSNP)) {
                                    // check whether the coding is identical
                                    if (allele1SNP.equals(allele1FirstSNP) && allele2SNP.equals(allele2FirstSNP)) {
                                        flipAlleles = false;
                                    } else if (allele1SNP.equals(allele2FirstSNP) && allele2SNP.equals(allele1FirstSNP)) {
                                        flipAlleles = true;
                                    } else {
                                        incompatible = true;
                                        System.out.println("Incompatible alleles for dataset: " + settings[d].name + "\t" + snp + "\t" + allele1FirstSNP + "/" + allele2FirstSNP +"("+minorAlleleFirstSNP+")\t" + allele1SNP + "/" + allele2SNP +"("+minorAlleleSNP+")");
                                    }
                                }
                            }

                            if (allele1SNP.equals(allele1FirstSNP) && allele2SNP.equals(allele2FirstSNP)) {
                                flipAlleles = false;
                            } else if (allele1SNP.equals(allele2FirstSNP) && allele2SNP.equals(allele1FirstSNP)) {
                                flipAlleles = true;
                            } else {
                                allele1SNP = BaseAnnot.getComplement(allele1SNP);
                                allele2SNP = BaseAnnot.getComplement(allele2SNP);

                                if (allele1SNP.equals(allele1FirstSNP) && allele2SNP.equals(allele2FirstSNP)) {
                                    flipAlleles = false;
                                } else if (allele1SNP.equals(allele2FirstSNP) && allele2SNP.equals(allele1FirstSNP)) {
                                    flipAlleles = true;
                                } else {
                                    System.out.println("Incompatible alleles for dataset: " + settings[d].name + "\t" + snp + "\t" + allele1FirstSNP + "/" + allele2FirstSNP + "\t" + allele1SNP + "/" + allele2SNP);
                                    incompatible = true;
                                }
                            }

                            if (!incompatible) {
                                if (!flipAlleles) {
                                    genotypeLabels[d][0] = BaseAnnot.toString(snpObj.getAlleles()[0]) + BaseAnnot.toString(snpObj.getAlleles()[0]);
                                    genotypeLabels[d][1] = BaseAnnot.toString(snpObj.getAlleles()[0]) + BaseAnnot.toString(snpObj.getAlleles()[1]);
                                    genotypeLabels[d][2] = BaseAnnot.toString(snpObj.getAlleles()[1]) + BaseAnnot.toString(snpObj.getAlleles()[1]);
                                } else {
                                    genotypeLabels[d][2] = BaseAnnot.toString(snpObj.getAlleles()[0]) + BaseAnnot.toString(snpObj.getAlleles()[0]);
                                    genotypeLabels[d][1] = BaseAnnot.toString(snpObj.getAlleles()[0]) + BaseAnnot.toString(snpObj.getAlleles()[1]);
                                    genotypeLabels[d][0] = BaseAnnot.toString(snpObj.getAlleles()[1]) + BaseAnnot.toString(snpObj.getAlleles()[1]);
                                }
                            }
                        }

                        if (!incompatible) {
                            String[] individuals = gt.getIndividuals();
                            int gt0ctr = 0;
                            int gt1ctr = 0;
                            int gt2ctr = 0;

                            for (int i = 0; i < individuals.length; i++) {
                                String ind = individuals[i];
                                if (gt.getIsIncluded()[i]) {
                                    String expSample = ds[d].getGenotypeToExpressionCouplings().get(ind);

                                    if (expSample != null) {
                                        int sampleId = ds[d].getExpressionData().getIndividualId(expSample);
                                        int genotype = snpObj.getGenotypes()[i];
                                        if (flipAlleles) {
                                            genotype = Math.abs(genotype - 2);
                                        }


                                        switch (genotype) {
                                            case 0:
                                                gt0ctr++;
                                                valsgt0.add(ds[d].getExpressionData().getMatrix()[probeId][sampleId]);
                                                break;
                                            case 1:
                                                gt1ctr++;
                                                valsgt1.add(ds[d].getExpressionData().getMatrix()[probeId][sampleId]);
                                                break;
                                            case 2:
                                                gt2ctr++;
                                                valsgt2.add(ds[d].getExpressionData().getMatrix()[probeId][sampleId]);
                                                break;
                                        }
                                    }
                                }
                            }


                            genotypeLabels[d][0] += " (" + gt0ctr + ")";
                            genotypeLabels[d][1] += " (" + gt1ctr + ")";
                            genotypeLabels[d][2] += " (" + gt2ctr + ")";


                            // put them in the array
                            vals[d][0] = copyToPrimitiveArr(valsgt0);
                            vals[d][1] = copyToPrimitiveArr(valsgt1);
                            vals[d][2] = copyToPrimitiveArr(valsgt2);


                        }

                        loader.close();

                    }

                }
                if (dsctr > 0) {

                    ViolinBoxPlot vbp = new ViolinBoxPlot();

                    // width: 3 genotypes * num datasets * 3 genotypes * marginbetween
                    int marginBetween = 5;
                    int margin = 10;

                    int width = (3 * ds.length) + (3 * marginBetween * ds.length) + 3 * 100;
                    int height = 3 * width;
                    System.out.println(width);
                    String outfilename = outdir + snp + "-" + probe + "-ViolinBoxPlot.pdf";
//                    vbp.draw(vals, datasetNames, genotypeLabels, width, height, margin, marginBetween, ViolinBoxPlot.Output.PDF, outfilename);
                }
            }

        }

    }

    private double[] copyToPrimitiveArr(ArrayList<Double> valsObj) {
        double[] vals = new double[valsObj.size()];
        for (int i = 0; i < valsObj.size(); i++) {
            vals[i] = valsObj.get(i);
        }
        return vals;
    }
}
