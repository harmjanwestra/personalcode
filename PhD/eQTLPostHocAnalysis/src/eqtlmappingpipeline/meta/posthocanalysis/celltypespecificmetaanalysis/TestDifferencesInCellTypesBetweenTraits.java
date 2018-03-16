/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import umcg.genetica.gwas.Independifier;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.gwas.Dependifier;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.TTest;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harm-jan
 */
public class TestDifferencesInCellTypesBetweenTraits {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            TestDifferencesInCellTypesBetweenTraits t = new TestDifferencesInCellTypesBetweenTraits();

            boolean eQTLsOnColumns = false;
            double gwaspvaluethreshold = 5E-8;
            double proxyldthreshold = 0.5;

            int maxproxydistance = 1000000;
            double pruneLDThreshold = 0.2;
            int minNumberOfTraitSNPs = 10;
            boolean filterEQTLsForGWASSNPs = false;
            boolean independifyInput = false;
//
//            String gwasCatalogLoc = "C:\\Work\\2013-05-22-gwascatalog.txt";
//            String matrixLoc = "C:\\Work\\CellTypeSpecificity\\table04.txt";
//            String proxyreferencedataset = "C:\\Work\\GGDs\\HapMapCEU\\TriTyper\\";
//            String outputdir = "C:\\Work\\CellTypeSpecificity\\";

            String gwasCatalogLoc = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/2013-05-22-gwascatalog.txt";
//            String matrixLoc = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-09-MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
            String matrixLoc = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/2013-07-11-DannyTVectorToInteractionTermCorrelationsPerCellType/table04.txt";
            // 
            String proxyreferencedataset = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap2r24-CEU/";
            String outputdir = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/2013-07-11-DannyTVectorToInteractionTermCorrelationsPerCellType/UnPruned";
            String metamatrix = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-09-MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
//            String metamatrix = null;
            t.run(gwasCatalogLoc, matrixLoc, eQTLsOnColumns, gwaspvaluethreshold, proxyldthreshold, maxproxydistance, proxyreferencedataset, pruneLDThreshold, minNumberOfTraitSNPs, outputdir, metamatrix, independifyInput, filterEQTLsForGWASSNPs);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String gwasCatalogLoc, String matrixLoc, boolean eQTLsOnColumns, double gwaspvaluethreshold,
            double proxyldthreshold, int maxproxydistance, String proxyreferencedataset, double pruneLDThreshold,
            int minNumberOfTraitSNPs, String outputdir, String metaMatrixfile, boolean independiyInput, boolean filterEQTLsForGWASSNPs) throws IOException {


        HashMap<String, Boolean> flipFx = null;
        if (metaMatrixfile != null) {
            flipFx = flipEffect(metaMatrixfile);
        }

        // load genotype datasets
        TriTyperGenotypeData genotypeData = new TriTyperGenotypeData(proxyreferencedataset);
        SNPLoader loader = genotypeData.createSNPLoader();

        // load GWAS catalog
        GWASCatalog c = new GWASCatalog(gwasCatalogLoc);
        GWASTrait[] traits = c.getTraits();
        Independifier independifier = new Independifier(genotypeData, loader);
        Dependifier proxyfinder = new Dependifier(genotypeData, loader);



        // procedure:
        // 1. load eqtl matrix
        // 2. iterate traits
        // 3. get trait snps per trait
        // 4. get trait proxies
        // 4. get the list of eQTLs for trait
        // 5. prune eQTL list for trait
        // 6. if eQTL list > 10 eQTLs, prune the rest of the eQTLs
        // 7. perform T-test on eQTL scores.

        // 1. load matrix
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(matrixLoc);
        if (!eQTLsOnColumns) {
            ds.transposeDataset();
        }

        List<String> allEQTLs = ds.colObjects;

        boolean[] excludeEQTL = new boolean[allEQTLs.size()];

        String[] alleQTLSNPsInDataset = new String[allEQTLs.size()];
        for (int i = 0; i < alleQTLSNPsInDataset.length; i++) {
            alleQTLSNPsInDataset[i] = allEQTLs.get(i).split("-")[0];
        }


        if (filterEQTLsForGWASSNPs) {
            System.out.println("Filtering eQTLs for GWAS SNPs..");
            // get all proxies for all GWAS SNPs..
            GWASSNP[] gwasSNPs = c.getSnpsArray();
            HashSet<String> allGWASSNPs = new HashSet<String>();
            for (int i = 0; i < gwasSNPs.length; i++) {
                GWASTrait[] traitsForSNP = gwasSNPs[i].getAssociatedTraitsArray();

                for (GWASTrait trait : traitsForSNP) {
                    Double p = gwasSNPs[i].getPValueAssociatedWithTrait(trait);
                    if (p != null && p < gwaspvaluethreshold) {
                        allGWASSNPs.add(gwasSNPs[i].getName());
                    }
                }

            }

            System.out.println("Detected: " + allGWASSNPs.size() + " SNPs with p < " + gwaspvaluethreshold);

            HashSet<String> allEQTLSNPsHash = new HashSet<String>();
            allEQTLSNPsHash.addAll(Arrays.asList(allEQTLSNPsHash.toArray(new String[0])));
            int ctr = 0;
            for (String gwaSSNP : allGWASSNPs) {
                if (allEQTLSNPsHash.contains(gwaSSNP)) {
                    ctr++;
                }
            }
            System.out.println("EQTL SNPs that are also GWAS SNPS: " + ctr);

            System.out.println("Checking whether any GWAS SNP is linked to any eQTL SNP.");
            HashSet<String> proxiesForGWASSNPs = proxyfinder.dependify(allGWASSNPs.toArray(new String[0]), alleQTLSNPsInDataset, proxyldthreshold, maxproxydistance);
            System.out.println("Found: " + proxiesForGWASSNPs.size() + " proxies linked to the eQTL SNPs.");


            System.out.println("Filtering the list of EQTLs. Removing eQTLs where the SNP is not a (proxy of a) GWAS SNP.");
            // now filter the eQTLs for GWAS SNPs..
            ctr = 0;
            for (int i = 0; i < alleQTLSNPsInDataset.length; i++) {
                String snp = alleQTLSNPsInDataset[i];
                if (!proxiesForGWASSNPs.contains(snp) && !allGWASSNPs.contains(snp)) {
                    excludeEQTL[i] = true;
                    ctr++;
                }
            }

            System.out.println(ctr + " / " + alleQTLSNPsInDataset.length + " eQTLs excluded after filtering for GWAS SNPs. Remaining eQTLs: " + (alleQTLSNPsInDataset.length - ctr));

        }

        System.out.println("Loaded data.. now processing " + traits.length + " traits");

        TextFile outputfile = new TextFile(outputdir + "PValuesPerTraitAndTissue.txt", TextFile.W);
        String header = "Pvalue\tPWilcoxon\tWilcoxonAUC\tTrait\tCellType\tTraitEQTLSNPs\tOtherEQTLSNPs";
        outputfile.writeln(header);
        // 2. iterate traits
        int ctr = 0;

        WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();

        for (GWASTrait trait : traits) {
            ctr++;
            GWASSNP[] traitSNPs = trait.getSNPs(gwaspvaluethreshold);

            // 3. get snps per trait
            HashSet<String> traitSNPStr = new HashSet<String>();
            for (GWASSNP snp : traitSNPs) {
                traitSNPStr.add(snp.getName());
            }

            // 4. get trait proxies
            HashSet<String> proxies = proxyfinder.dependify(traitSNPStr.toArray(new String[0]), proxyldthreshold, maxproxydistance);

            HashSet<String> traitEQTLSNPs = new HashSet<String>();
            HashSet<String> otherEQTLSNPs = new HashSet<String>();

            for (int i = 0; i < allEQTLs.size(); i++) {
                String snp = alleQTLSNPsInDataset[i];
                if (!excludeEQTL[i]) {
                    if (proxies.contains(snp)) {
                        traitEQTLSNPs.add(snp);
                    } else {
                        otherEQTLSNPs.add(snp);
                    }
                }
            }

            System.out.println(ctr + "/" + traits.length + "\t" + trait.getName() + "\tNrSNPs: " + traitSNPs.length + "\tProxies: " + proxies.size() + "\tNrEQTLSNPs: " + traitEQTLSNPs.size() + "\tRemainingSNPs: " + otherEQTLSNPs.size());

            if (traitEQTLSNPs.size() >= minNumberOfTraitSNPs) {
                String[] independifiedTraitEQTLSNPs = traitEQTLSNPs.toArray(new String[0]);
                String[] independifiedOtherEQTLSNPs = otherEQTLSNPs.toArray(new String[0]);
                if (independiyInput) {
                    independifiedTraitEQTLSNPs = independifier.independify(independifiedTraitEQTLSNPs, pruneLDThreshold, maxproxydistance);
                    independifiedOtherEQTLSNPs = independifier.independify(independifiedOtherEQTLSNPs, pruneLDThreshold, maxproxydistance);
                }


                if (independifiedTraitEQTLSNPs.length >= minNumberOfTraitSNPs) {
                    HashSet<String> independifiedTraitEQTLSNPsHash = new HashSet<String>();
                    HashSet<String> independifiedOtherEQTLSNPsHash = new HashSet<String>();
                    independifiedTraitEQTLSNPsHash.addAll(Arrays.asList(independifiedTraitEQTLSNPs));
                    independifiedOtherEQTLSNPsHash.addAll(Arrays.asList(independifiedOtherEQTLSNPs));

                    for (int row = 0; row < ds.nrRows; row++) {
                        String rowName = ds.rowObjects.get(row);
                        ArrayList<Double> traitEQTLSNPVals = new ArrayList<Double>();
                        ArrayList<Double> otherEQTLSNPVals = new ArrayList<Double>();
                        for (int col = 0; col < ds.nrCols; col++) {
                            if (!excludeEQTL[col]) {
                                String eQTL = allEQTLs.get(col);
                                String snp = alleQTLSNPsInDataset[col];

                                double val = ds.rawData[row][col];


                                if (flipFx != null && !rowName.equals("MainEffectZScore")) {
                                    boolean flip = flipFx.get(eQTL);
                                    if (flip) {
                                        val *= -1;
                                    }
                                }

                                if (independifiedTraitEQTLSNPsHash.contains(snp)) {
                                    traitEQTLSNPVals.add(val);
                                }
                                if (independifiedOtherEQTLSNPsHash.contains(snp)) {
                                    otherEQTLSNPVals.add(val);
                                }
                            }

                        }
                        double[] xarr = Primitives.toPrimitiveArr(traitEQTLSNPVals.toArray(new Double[0]));
                        double[] yarr = Primitives.toPrimitiveArr(otherEQTLSNPVals.toArray(new Double[0]));
                        double p = TTest.test(xarr, yarr);

                        double pwilcoxon = mwm.returnWilcoxonMannWhitneyPValue(xarr, yarr);
                        double auc = mwm.getAUC();

                        String output = p + "\t" + pwilcoxon + "\t" + auc + "\t" + trait.getName() + "\t" + ds.rowObjects.get(row) + "\t" + traitEQTLSNPVals.size() + "\t" + otherEQTLSNPVals.size();
                        System.out.println(output);
                        outputfile.writeln(output);

                    }

                }
            }
        }
        outputfile.close();
        // close the genotype loader
        loader.close();
    }

    private HashMap<String, Boolean> flipEffect(String metaMatrixLoc) throws IOException {
        DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(metaMatrixLoc);
        HashMap<String, Boolean> output = new HashMap<String, Boolean>();

        Integer rowIdForInteractionTerm = metaMatrix.hashRows.get("CellTypeInteractionZScore");
        Integer rowIdForMainEffect = metaMatrix.hashRows.get("MainEffectZScore");

        for (int col = 0; col < metaMatrix.nrCols; col++) {

            double main = metaMatrix.rawData[rowIdForMainEffect][col];
            double inte = metaMatrix.rawData[rowIdForInteractionTerm][col];


            if (main < 0) {
                output.put(metaMatrix.colObjects.get(col), true);
            } else {
                output.put(metaMatrix.colObjects.get(col), false);
            }

        }

        return output;
    }
}
