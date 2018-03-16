/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.gwas.Independifier;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.gwas.Dependifier;
import umcg.genetica.io.Gpio;
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

            boolean eQTLsOnColumns = true;
            double gwaspvaluethreshold = 5E-8;
            double proxyldthreshold = 2;

            int maxproxydistance = 1000000;
            double pruneLDThreshold = 0.2;
            int minNumberOfTraitSNPs = 10;
            boolean filterEQTLsForGWASSNPs = false;
            boolean independifyInput = true;
            boolean flipEffectsAccordingToMetaMatrixMainEffectDirection = true;
//
//            String gwasCatalogLoc = "C:\\Work\\2013-05-22-gwascatalog.txt";
//            String matrixLoc = "C:\\Work\\CellTypeSpecificity\\table04.txt";
//            String proxyreferencedataset = "C:\\Work\\GGDs\\HapMapCEU\\TriTyper\\";
//            String outputdir = "C:\\Work\\CellTypeSpecificity\\";


            // Meta w/ SHIP: SHIP+Rotterdam calculated by cohorts themselves
//            String gwasCatalogLoc = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-07-29-hjcatalog.txt";
            String gwasCatalogLoc = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2013-07-22-gwascatalog.txt";
            String proxyreferencedataset = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap2r24-CEU/";

//            for (int gwasf = 0; gwasf < 2; gwasf++) {
//                String gwasfStr = "AllSNPs";
//                filterEQTLsForGWASSNPs = false;
//                if (gwasf == 1) {
//                    filterEQTLsForGWASSNPs = true;
//                    gwasfStr = "GWASSNPs";
//                }
//
//                for (int prune = 0; prune < 2; prune++) {
//
//                    String pruneStr = "Unpruned";
//                    independifyInput = false;
//                    if (prune == 1) {
//                        independifyInput = true;
//                        pruneStr = "Pruned";
//                    }
//
//                    for (int dataset = 0; dataset < 2; dataset++) {
//                        String matrixLoc = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-24-12KEQTLs-WSHIP-TREND/MetaAnalysis/MetaAnalysisZScoreMatrix-CorrelatedWithCellCountCorrelations.txt";
//                        String metamatrix = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-24-12KEQTLs-WSHIP-TREND/MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
//                        String outputdir = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-29-FakeGWAS/MetaWSHIP-TREND/InteractionVectorsCorrelatedWithCellCounts/" + pruneStr + "-" + gwasfStr + "-TTestsForTraitsUsingCellCountCorrelations";
//
//                        if (dataset == 1) {
//                            matrixLoc = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/MetaAnalysis/CellCountCorrelationMatrixCorrelatedWithInteractionTerms.txt";
//                            metamatrix = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
//                            outputdir = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-29-FakeGWAS/MetaWOSHIP-TREND/" + pruneStr + "-" + gwasfStr + "-TTestsForTraitsUsingCellCountCorrelations";
//                        }
//
//                        Gpio.createDir(outputdir);
//                        t.testDifferenceWithinCellType(outputdir, gwaspvaluethreshold, proxyreferencedataset, proxyldthreshold, pruneLDThreshold, maxproxydistance, minNumberOfTraitSNPs, independifyInput, flipEffectsAccordingToMetaMatrixMainEffectDirection, metamatrix, eQTLsOnColumns, matrixLoc, gwasCatalogLoc, filterEQTLsForGWASSNPs);
////                        t.testDifferenceBetweenCellTypes(outputdir, gwaspvaluethreshold, proxyreferencedataset, proxyldthreshold, pruneLDThreshold, maxproxydistance, minNumberOfTraitSNPs, independifyInput, flipEffectsAccordingToMetaMatrixMainEffectDirection, metamatrix, eQTLsOnColumns, matrixLoc, gwasCatalogLoc, filterEQTLsForGWASSNPs);
//                    }
//                }
//            }

            String matrixLoc = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-31-MetaAnalysis5Datasets/MetaAnalysis/CellTypeSpecificityMatrix.binary";
            String metamatrix = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-24-12KEQTLs-WSHIP-TREND/MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
            String outputdir = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-31-MetaAnalysis5Datasets/TTestsForTraitsInteractionZScores/TTestsForTraitsUsingCellCountCorrelationsTop250";
            Gpio.createDir(outputdir);
            t.testDifferenceWithinCellType(outputdir, gwaspvaluethreshold, proxyreferencedataset, proxyldthreshold, pruneLDThreshold, maxproxydistance, minNumberOfTraitSNPs, independifyInput, flipEffectsAccordingToMetaMatrixMainEffectDirection, metamatrix, eQTLsOnColumns, matrixLoc, gwasCatalogLoc, filterEQTLsForGWASSNPs);


            System.exit(0);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    private HashMap<String, Double> mainEffectPerEQTL;
    private GWASCatalog c;
    private GWASTrait[] traits;
    private Independifier independifier;
    private Dependifier proxyfinder;
    private DoubleMatrixDataset<String, String> ds;
    private List<String> allEQTLs;
    private HashSet<String> eQTLsToInclude;
    private String[] alleQTLSNPsInDataset;
    private HashMap<String, Boolean> flipFx;
    private WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();
    private String outputdir;

    public void init(boolean flipEffectsAccordingToMetaMatrixMainEffectDirection, String metaMatrixfile, boolean eQTLsOnColumns, String matrixLoc, String gwasCatalogLoc,
            TriTyperGenotypeData genotypeData, SNPLoader loader, boolean filterEQTLsForGWASSNPs, double gwaspvaluethreshold) throws IOException {
        flipFx = null;

        if (flipEffectsAccordingToMetaMatrixMainEffectDirection) {
            flipFx = flipEffect(metaMatrixfile);
        }
        mainEffectPerEQTL = getMainEffectZScorePerEQTL(metaMatrixfile);

        // load genotype datasets


        // load GWAS catalog
        c = new GWASCatalog(gwasCatalogLoc);
        traits = c.getTraits();

        independifier = new Independifier(genotypeData, loader);
        proxyfinder = new Dependifier(genotypeData, loader);


        // 1. load matrix
        HashSet<String> set = new HashSet<String>();
//        set.add("CellTypeInteractionZScore");
        set.add("6270612");set.add("3990112");set.add("5810435");set.add("2030484");set.add("5690187");set.add("2120152");set.add("7550292");set.add("7040035");set.add("3140750");set.add("3780452");set.add("4280603");set.add("1820528");set.add("20091");set.add("6650594");set.add("1070242");set.add("780402");set.add("3850561");set.add("6520451");set.add("2710653");set.add("7330671");set.add("3840400");set.add("830164");set.add("1780356");set.add("5050681");set.add("1850041");set.add("3420372");set.add("6180048");set.add("5900468");set.add("5080154");set.add("7160202");set.add("6760291");set.add("360187");set.add("1090390");set.add("4780048");set.add("1980768");set.add("5870594");set.add("1780709");set.add("4900168");set.add("CellTypeZScore");set.add("940386");set.add("10142");set.add("7200603");set.add("1030348");set.add("6420671");set.add("6650161");set.add("2260551");set.add("1710278");set.add("6110091");set.add("3290091");set.add("1090692");set.add("7200315");set.add("5900129");set.add("1990110");set.add("6100767");set.add("CellTypeInteractionZScore");set.add("5220189");set.add("50192");set.add("2480192");set.add("5890138");set.add("6110747");set.add("7210280");set.add("5270148");set.add("4900333");set.add("520706");set.add("620450");set.add("3460278");set.add("1240446");set.add("2600204");set.add("6550703");set.add("6660286");set.add("6280008");set.add("430356");set.add("2760259");set.add("1050482");set.add("1850220");set.add("6250110");set.add("2120689");set.add("1050215");set.add("5360273");set.add("1990717");set.add("7380626");set.add("5560400");set.add("1470196");set.add("2690541");set.add("160494");set.add("1510139");set.add("460333");set.add("3460201");set.add("3170482");set.add("1030431");set.add("6380717");set.add("5340338");set.add("3180465");set.add("3420719");set.add("6420242");set.add("3940026");set.add("6550600");set.add("6660674");set.add("6900079");set.add("6860347");set.add("4230097");set.add("5860215");set.add("5670465");set.add("3190634");set.add("540091");set.add("670086");set.add("160026");set.add("650634");set.add("4670327");set.add("2480468");set.add("2760255");set.add("5260288");set.add("990735");set.add("3120075");set.add("6860196");set.add("6110386");set.add("6370678");set.add("1820008");set.add("7320246");set.add("2760427");set.add("2690239");set.add("5270138");set.add("1070477");set.add("1690240");set.add("6330487");set.add("3290369");set.add("1010521");set.add("540390");set.add("3360634");set.add("6280458");set.add("1010253");set.add("6370435");set.add("6520605");set.add("6840468");set.add("7550112");set.add("4560730");set.add("6040008");set.add("2350703");set.add("5670735");set.add("6960274");set.add("3400482");set.add("5420538");set.add("620047");set.add("2570131");set.add("7650053");set.add("6180544");set.add("5420564");set.add("2640274");set.add("6940634");set.add("6200408");set.add("6590441");set.add("4880632");set.add("3370768");set.add("6040468");set.add("5820528");set.add("7570367");set.add("2600634");set.add("6520646");set.add("1090221");set.add("3060692");set.add("5360719");set.add("830093");set.add("3840437");set.add("130392");set.add("6370113");set.add("4670102");set.add("3780100");set.add("2340500");set.add("7200242");set.add("670273");set.add("2630687");set.add("6100424");set.add("10594");set.add("4210129");set.add("3930343");set.add("6420446");set.add("2850184");set.add("5870632");set.add("2360630");set.add("5910113");set.add("1770730");set.add("3610768");set.add("160537");set.add("6200577");set.add("4260253");set.add("4850168");set.add("2030767");set.add("7510424");set.add("2060286");set.add("1990468");set.add("1780673");set.add("2450647");set.add("6860048");set.add("5890661");set.add("2340347");set.add("4490242");set.add("7320021");set.add("1030376");set.add("110255");set.add("3400332");set.add("2480170");set.add("1820324");set.add("4390022");set.add("3400392");set.add("6020575");set.add("6960209");set.add("4860563");set.add("20288");set.add("1240064");set.add("4010564");set.add("4150543");set.add("1940300");set.add("6660435");set.add("5900497");set.add("2140020");set.add("4060494");set.add("CellTypeSNPZScore");set.add("5220240");set.add("1110300");set.add("6580164");set.add("5050025");set.add("1410221");set.add("10685");set.add("1030209");set.add("5560020");set.add("110064");set.add("6100703");set.add("7000377");set.add("20142");set.add("460373");set.add("3390048");set.add("670470");set.add("6510601");set.add("10379");set.add("2450082");set.add("770561");set.add("1260136");set.add("1450608");set.add("7210450");set.add("4180717");set.add("1190327");set.add("6020681");set.add("2190411");set.add("5130750");set.add("7200358");set.add("70403");set.add("3060397");set.add("4230373");set.add("2030068");set.add("7040411");
        ds = new DoubleMatrixDataset<String, String>(matrixLoc, set);
        if (!eQTLsOnColumns) {
            ds.transposeDataset();
        }

        ds.save(outputdir + "LoadedMatrix.txt");

        allEQTLs = ds.colObjects;
        eQTLsToInclude = null;
        alleQTLSNPsInDataset = new String[allEQTLs.size()];
        for (int i = 0; i < alleQTLSNPsInDataset.length; i++) {
            String eQTL = allEQTLs.get(i);
            String[] eQTLElems = eQTL.split("-");
            alleQTLSNPsInDataset[i] = eQTLElems[0];
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
            System.out.println("Filtering the list of EQTLs. Removing eQTLs where the SNP is not a (proxy of a) GWAS SNP.");
            // now filter the eQTLs for GWAS SNPs..
            int ctr = 0;
            eQTLsToInclude = new HashSet<String>();
            for (int i = 0; i < allEQTLs.size(); i++) {
                String snp = allEQTLs.get(i).split("-")[0];
                if (allGWASSNPs.contains(snp)) {
                    eQTLsToInclude.add(allEQTLs.get(i));
                    ctr++;
                }
            }
            System.out.println(ctr + " / " + alleQTLSNPsInDataset.length + " eQTLs included after filtering for GWAS SNPs.");
        }
    }

    public void testDifferenceWithinCellType(String outputdir, double gwaspvaluethreshold, String referenceDatasetLocation, double proxyldthreshold, double pruneLDThreshold, int maxproxydistance,
            int minNumberOfTraitSNPs, boolean independiyInput,
            boolean flipEffectsAccordingToMetaMatrixMainEffectDirection, String metaMatrixFile, boolean eQTLsOnColumns, String matrixLoc, String gwasCatalogLoc, boolean filtereQTLsForGWASSNPs) throws IOException {

        TriTyperGenotypeData genotypeData = new TriTyperGenotypeData(referenceDatasetLocation);
        SNPLoader loader = genotypeData.createSNPLoader();
        this.outputdir = outputdir;

        init(flipEffectsAccordingToMetaMatrixMainEffectDirection, metaMatrixFile, eQTLsOnColumns, matrixLoc, gwasCatalogLoc, genotypeData, loader, filtereQTLsForGWASSNPs, gwaspvaluethreshold);




        String fileSuffix = "";
        String fileName1 = outputdir + "PValuesPerTraitAndTissue.txt";
        String fileName2 = outputdir + "PValuesPerTraitAndTissueCorrelationsPerTraitEQTL.txt";

//        if (Gpio.exists(fileName1)) {
//            int suffixCtr = 0;
//            while (Gpio.exists(fileName1)) {
//                fileName1 = outputdir + "PValuesPerTraitAndTissue." + suffixCtr + ".txt";
//                fileName2 = outputdir + "PValuesPerTraitAndTissueCorrelationsPerTraitEQTL." + suffixCtr + ".txt";
//            }
//        }
        TextFile outputfile = new TextFile(fileName1, TextFile.W);
        String header = "Pvalue\tPWilcoxon\tWilcoxonAUC\tTrait\tCellType\tTraitEQTLSNPs\tOtherEQTLSNPs";

        TextFile outputFileCorelations = new TextFile(fileName2, TextFile.W);
        outputFileCorelations.writeln("Trait\teQTL\tCellCount\tCorrelation");
        outputfile.writeln(header);
        // 2. iterate traits
        int ctr = 0;
        System.out.println("Loaded data.. now processing " + traits.length + " traits");


        for (GWASTrait trait : traits) {
            ctr++;
            GWASSNP[] traitSNPs = trait.getSNPs(gwaspvaluethreshold);

            // 3. get snps per trait
            HashSet<String> traitSNPStrHash = new HashSet<String>();
            for (GWASSNP snp : traitSNPs) {
                traitSNPStrHash.add(snp.getName());
            }

            // 4. get trait proxies
            HashSet<String> proxies = proxyfinder.dependify(traitSNPStrHash.toArray(new String[0]), proxyldthreshold, maxproxydistance);
            traitSNPStrHash.addAll(proxies);

            HashSet<String> traitEQTLSNPs = new HashSet<String>();
            HashSet<String> otherEQTLSNPs = new HashSet<String>();

            ArrayList<String> allTraitEQTLs = new ArrayList<String>();
            for (int i = 0; i < allEQTLs.size(); i++) {
                String snp = alleQTLSNPsInDataset[i];
                if (eQTLsToInclude == null || eQTLsToInclude.contains(allEQTLs.get(i))) {
                    if (traitSNPStrHash.contains(snp)) {
                        allTraitEQTLs.add(allEQTLs.get(i));
                        traitEQTLSNPs.add(snp);
                    } else {
                        otherEQTLSNPs.add(snp);
                    }
                }
            }

            System.out.println(ctr + "/" + traits.length + "\t" + trait.getName() + "\tNrSNPs: " + traitSNPs.length + "\tProxies: " + proxies.size() + "\tNrEQTLSNPs: " + traitEQTLSNPs.size() + "\tNr Trait eQTLs: " + allTraitEQTLs.size() + "\tRemainingSNPs: " + otherEQTLSNPs.size() + "\tRemaining eQTLs: " + (allEQTLs.size() - allTraitEQTLs.size()));

            if (traitEQTLSNPs.size() >= minNumberOfTraitSNPs) {
                String[] independifiedTraitEQTLSNPs = traitEQTLSNPs.toArray(new String[0]);
                String[] independifiedOtherEQTLSNPs = otherEQTLSNPs.toArray(new String[0]);
                if (independiyInput) {
                    // independification returns a ; separated list of SNPs.
                    independifiedTraitEQTLSNPs = independifier.independify(independifiedTraitEQTLSNPs, pruneLDThreshold, maxproxydistance);
                    //independifiedOtherEQTLSNPs = independifier.independify(independifiedOtherEQTLSNPs, pruneLDThreshold, maxproxydistance);
                }

                if (independifiedTraitEQTLSNPs.length >= minNumberOfTraitSNPs) {

                    // split the independified trait-eQTL SNPs, take the ones with the highest meta-analysis Z-score
                    HashSet<String> independifiedTraitEQTLSNPsHash = selectTopSNPsFromIndependifiedClump(independifiedTraitEQTLSNPs, allEQTLs, mainEffectPerEQTL);
                    HashSet<String> independifiedOtherEQTLSNPsHash = selectTopSNPsFromIndependifiedClump(independifiedOtherEQTLSNPs, allEQTLs, mainEffectPerEQTL);

                    double[][] allYVals1 = new double[ds.nrRows][0];
                    double[][] allYVals2 = new double[ds.nrRows][0];
                    int fullYLength1 = 0;
                    int fullYLength2 = 0;
                    for (int row = 0; row < ds.nrRows; row++) {
                        String cellCountName = ds.rowObjects.get(row);
                        ArrayList<Double> traitEQTLSNPVals = new ArrayList<Double>();
                        ArrayList<Double> otherEQTLSNPVals = new ArrayList<Double>();
                        for (int col = 0; col < ds.nrCols; col++) {
                            String eQTL = allEQTLs.get(col);
                            if (eQTLsToInclude == null || eQTLsToInclude.contains(eQTL)) {
                                String snp = alleQTLSNPsInDataset[col];
                                double val = ds.rawData[row][col];
                                if (flipFx != null && !cellCountName.equals("MainEffectZScore")) {
                                    boolean flip = flipFx.get(eQTL);
                                    if (flip) {
                                        val *= -1;
                                    }
                                }

                                if (independifiedTraitEQTLSNPsHash.contains(snp)) {
                                    traitEQTLSNPVals.add(val);
                                    outputFileCorelations.writeln(trait.getName() + "\t" + eQTL + "\t" + ds.rowObjects.get(row) + "\t" + val);
                                }
                                if (independifiedOtherEQTLSNPsHash.contains(snp)) {
                                    otherEQTLSNPVals.add(val);
                                }
                            }

                        }
                        double[] xarr = Primitives.toPrimitiveArr(traitEQTLSNPVals.toArray(new Double[0]));

                        double[] yarr = Primitives.toPrimitiveArr(otherEQTLSNPVals.toArray(new Double[0]));
                        allYVals1[row] = xarr;
                        fullYLength1 += xarr.length;
                        allYVals2[row] = yarr;
                        fullYLength2 += yarr.length;
                        double p = TTest.test(xarr, yarr);

                        double pwilcoxon = mwm.returnWilcoxonMannWhitneyPValue(xarr, yarr);
                        double auc = mwm.getAUC();

                        String output = p + "\t" + pwilcoxon + "\t" + auc + "\t" + trait.getName() + "\t" + ds.rowObjects.get(row) + "\t" + traitEQTLSNPVals.size() + "\t" + otherEQTLSNPVals.size();
                        System.out.println(output);
                        outputfile.writeln(output);

                    }
                    double[] plotY = new double[fullYLength1 + fullYLength2];
                    double[] plotX = new double[fullYLength1 + fullYLength2];
                    int ctr1 = 0;
                    for (int i = 0; i < allYVals1.length; i++) {
                        for (int j = 0; j < allYVals1[i].length; j++) {
                            plotX[ctr1] = (i + 0.25) + (Math.random() * 0.125);
                            plotY[ctr1] = allYVals1[i][j];
                            ctr1++;
                        }
                        for (int j = 0; j < allYVals2[i].length; j++) {
                            plotX[ctr1] = (i + 0.75) + (Math.random() * 0.125);
                            plotY[ctr1] = allYVals2[i][j];
                            ctr1++;
                        }
                    }
                    new ScatterPlot(1500, 500, plotX, plotY, ScatterPlot.OUTPUTFORMAT.PNG, outputdir + trait.getName() + ".png");

                }
            }
        }
        outputFileCorelations.close();
        outputfile.close();
        // close the genotype loader
        loader.close();
    }

    private void testDifferenceBetweenCellTypes(String outputdir, double gwaspvaluethreshold, String referenceDatasetLocation, double proxyldthreshold, double pruneLDThreshold, int maxproxydistance,
            int minNumberOfTraitSNPs, boolean independiyInput,
            boolean flipEffectsAccordingToMetaMatrixMainEffectDirection, String metaMatrixFile, boolean eQTLsOnColumns, String matrixLoc, String gwasCatalogLoc, boolean filtereQTLsForGWASSNPs) throws IOException {
        TriTyperGenotypeData genotypeData = new TriTyperGenotypeData(referenceDatasetLocation);
        SNPLoader loader = genotypeData.createSNPLoader();

        init(flipEffectsAccordingToMetaMatrixMainEffectDirection, metaMatrixFile, eQTLsOnColumns, matrixLoc, gwasCatalogLoc, genotypeData, loader, filtereQTLsForGWASSNPs, gwaspvaluethreshold);

        System.out.println("Loaded data.. now processing " + traits.length + " traits");

        String fileSuffix = "";
        String fileName1 = outputdir + "PValuesTestingBetweenTraits.txt";
        String fileName2 = outputdir + "PValuesTestingBetweenTraitsCorrelationsPerTraitEQTL.txt";

        if (Gpio.exists(fileName1)) {
            int suffixCtr = 0;
            while (Gpio.exists(fileName1)) {
                fileName1 = outputdir + "PValuesTestingBetweenTraits." + suffixCtr + ".txt";
                fileName2 = outputdir + "PValuesTestingBetweenTraitsCorrelationsPerTraitEQTL." + suffixCtr + ".txt";
            }
        }
        TextFile outputfile = new TextFile(fileName1, TextFile.W);
        String header = "Pvalue\tPWilcoxon\tWilcoxonAUC\tTrait\tCellType\tTraitEQTLSNPs\tOtherEQTLSNPs";

        TextFile outputFileCorelations = new TextFile(fileName2, TextFile.W);
        outputFileCorelations.writeln("Trait\teQTL\tCellCount\tCorrelation");
        outputfile.writeln(header);
        // 2. iterate traits
        int ctr = 0;

        for (GWASTrait trait : traits) {
            ctr++;
            GWASSNP[] traitSNPs = trait.getSNPs(gwaspvaluethreshold);

            // 3. get snps per trait
            HashSet<String> traitSNPStrHash = new HashSet<String>();
            for (GWASSNP snp : traitSNPs) {
                traitSNPStrHash.add(snp.getName());
            }

            // 4. get trait proxies
            HashSet<String> proxies = proxyfinder.dependify(traitSNPStrHash.toArray(new String[0]), proxyldthreshold, maxproxydistance);
            traitSNPStrHash.addAll(proxies);

            HashSet<String> traitEQTLSNPs = new HashSet<String>();
//            HashSet<String> otherEQTLSNPs = new HashSet<String>();

            ArrayList<String> allTraitEQTLs = new ArrayList<String>();
            for (int i = 0; i < allEQTLs.size(); i++) {
                String snp = alleQTLSNPsInDataset[i];
                if (eQTLsToInclude == null || eQTLsToInclude.contains(allEQTLs.get(i))) {
                    if (traitSNPStrHash.contains(snp)) {
                        allTraitEQTLs.add(allEQTLs.get(i));
                        traitEQTLSNPs.add(snp);
                    } else {
//                        otherEQTLSNPs.add(snp);
                    }
                }
            }

            if (traitEQTLSNPs.size() >= minNumberOfTraitSNPs) {
                String[] independifiedTraitEQTLSNPs = traitEQTLSNPs.toArray(new String[0]);
//                String[] independifiedOtherEQTLSNPs = otherEQTLSNPs.toArray(new String[0]);
                if (independiyInput) {
                    // independification returns a ; separated list of SNPs.
                    independifiedTraitEQTLSNPs = independifier.independify(independifiedTraitEQTLSNPs, pruneLDThreshold, maxproxydistance);
//                    independifiedOtherEQTLSNPs = independifier.independify(independifiedOtherEQTLSNPs, pruneLDThreshold, maxproxydistance);
                }

                if (independifiedTraitEQTLSNPs.length >= minNumberOfTraitSNPs) {
                    HashSet<String> independifiedTraitEQTLSNPsHash = selectTopSNPsFromIndependifiedClump(independifiedTraitEQTLSNPs, allEQTLs, mainEffectPerEQTL);
                    for (int row = 0; row < ds.nrRows; row++) {
                        String rowName = ds.rowObjects.get(row); // this is where the cell count names reside..
                        ArrayList<Double> traitEQTLSNPVals = new ArrayList<Double>();
                        ArrayList<Double> otherEQTLSNPVals = new ArrayList<Double>();

                        for (int col = 0; col < ds.nrCols; col++) {
                            String eQTL = allEQTLs.get(col);
                            if (eQTLsToInclude == null || eQTLsToInclude.contains(eQTL)) {
                                String snp = alleQTLSNPsInDataset[col];
                                if (independifiedTraitEQTLSNPsHash.contains(snp)) {
                                    int flip = 1;
                                    if (flipFx != null && !rowName.equals("MainEffectZScore")) {
                                        boolean flipbool = flipFx.get(eQTL);
                                        if (flipbool) {
                                            flip = -1;
                                        }
                                        for (int row2 = 0; row2 < ds.nrRows; row2++) {
                                            double val = ds.rawData[row2][col];
                                            if (row2 != row) {
                                                otherEQTLSNPVals.add(flip * val);
                                            } else {
                                                traitEQTLSNPVals.add(flip * val);
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // test
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
        loader.close();
    }

    private HashMap<String, Boolean> flipEffect(String metaMatrixLoc) throws IOException {

        HashSet<String> set = new HashSet<String>();
        set.add("CellTypeInteractionZScore");
        set.add("MainEffectZScore");
        DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(metaMatrixLoc, set);


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

    private HashMap<String, Double> getMainEffectZScorePerEQTL(String metaMatrixLoc) throws IOException {
        HashSet<String> set = new HashSet<String>();
        set.add("MainEffectZScore");
        DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(metaMatrixLoc, set);
        HashMap<String, Double> output = new HashMap<String, Double>();

        Integer rowIdForMainEffect = metaMatrix.hashRows.get("MainEffectZScore");
        for (int col = 0; col < metaMatrix.nrCols; col++) {

            double main = metaMatrix.rawData[rowIdForMainEffect][col];
            output.put(metaMatrix.colObjects.get(col), main);
        }

        return output;
    }

    public HashSet<String> selectTopSNPsFromIndependifiedClump(String[] independifiedClumps, List<String> allEQTLs, HashMap<String, Double> mainEffectPerEQTL) {
        HashSet<String> output = new HashSet<String>();
        for (String clump : independifiedClumps) {
            String[] clumpElems = clump.split(";");

            if (clumpElems.length == 1) {
                // there is only one snp in this clump.
                output.add(clumpElems[0]);
            } else {
                HashSet<String> clumpedSNPs = new HashSet<String>(Arrays.asList(clumpElems));
                double maxZ = -1;
                String selectedSNP = null;
                for (int i = 0; i < allEQTLs.size(); i++) {
                    String eQTL = allEQTLs.get(i);
                    String snp = eQTL.split("-")[0];
                    if (clumpedSNPs.contains(snp)) {
                        Double Z = mainEffectPerEQTL.get(eQTL);
                        if (Z != null) {
                            Z = Math.abs(Z);
                            if (selectedSNP == null || Z > maxZ) {
                                selectedSNP = snp;
                                maxZ = Z;
                            }
                        } else {
                            System.err.println("ERROR: Z is null for " + eQTL);
                        }
                    }
                }
                output.add(selectedSNP);
            }
        }
        return output;
    }
}
