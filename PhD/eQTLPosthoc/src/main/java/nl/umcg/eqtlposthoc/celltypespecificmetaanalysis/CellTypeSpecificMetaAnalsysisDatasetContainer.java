/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CellTypeSpecificMetaAnalsysisDatasetContainer {

    private final HashSet<String> allSNPs;
    private final HashSet<String> allEQTL;
    private final ArrayList<DoubleMatrixDataset<String, String>> datasets;
    private final SNP[][] snps;
    private final Integer[][] eQTLIndex;
    private final HashMap<String, Integer> eQTLToId;
    private final Integer[][] snpIndex;
    private final HashMap<String, Integer> snpToId;
    private final String[] allCovariateArr;
    private final Integer[][] covariateIndex;
    private final HashSet<String> allCovariates;
    private final HashMap<String, Integer> covariateToId;

    public CellTypeSpecificMetaAnalsysisDatasetContainer(String[] datasetFileDirs, boolean[] isHT12v3, HashMap<String, String> ht12v4ToHT12v3) throws IOException {
        this(datasetFileDirs, isHT12v3, ht12v4ToHT12v3, null);
    }

    public CellTypeSpecificMetaAnalsysisDatasetContainer(String[] datasetFileDirs, boolean[] isHT12v3, HashMap<String, String> ht12v4ToHT12v3, HashSet<String> requiredRows) throws IOException {
        // load the datasets..
        snps = new SNP[datasetFileDirs.length][0]; // JAGGED ARRAY FOR SNP STORAGE

        datasets = new ArrayList<DoubleMatrixDataset<String, String>>();
        allSNPs = new HashSet<String>();
        allEQTL = new HashSet<String>();
        for (int dsN = 0; dsN < datasetFileDirs.length; dsN++) {
            String dir = datasetFileDirs[dsN];
            // load the dataset....
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(dir + "CellTypeSpecificityMatrix.binary", requiredRows);
            datasets.add(ds);

            List<String> eqtls = ds.colObjects;

            for (int i = 0; i < eqtls.size(); i++) {
                String eQTL = eqtls.get(i);
                if (isHT12v3[dsN]) {
                    allEQTL.add(eQTL);
                } else {
                    String[] eqtlelems = eQTL.split("-");
                    eqtlelems[1] = ht12v4ToHT12v3.get(eqtlelems[1]);
                    allEQTL.add(Strings.concat(eqtlelems, Strings.dash));
                }
            }

            TextFile tfIn = new TextFile(dir + "SNPSummaryStatistics.txt", TextFile.R);
            tfIn.readLine(); // skip the header..
            ArrayList<SNP> allSNPsInDs = new ArrayList<SNP>();
            String[] elems = tfIn.readLineElems(TextFile.tab);
            HashSet<String> visitedSNPs = new HashSet<String>();
            while (elems != null) {
                // SNP	Chr	ChrPos	Alleles	MinorAllele	MAF	CallRate	HWE	GenotypesCalled
                String snp = elems[0];
                if (!visitedSNPs.contains(snp)) {
                    String alleles = elems[3];
                    String alleleAssessed = elems[4];
                    String maf = elems[5];
                    String cr = elems[6];
                    String hwe = elems[7];
                    String nrCalled = elems[8];

                    SNP s = new SNP();
                    s.setName(snp);
                    byte[] allelesB = new byte[2];
                    String[] alleleElems = alleles.split("/");
                    allelesB[0] = BaseAnnot.toByte(alleleElems[0]);
                    allelesB[1] = BaseAnnot.toByte(alleleElems[1]);
                    s.setAlleleCodes(allelesB);
                    byte alleleAssessedB = BaseAnnot.toByte(alleleAssessed);
                    s.setMinorAllele(alleleAssessedB);

                    try {
                        s.setMAF(Double.parseDouble(maf));
                        s.setCR(Double.parseDouble(cr));
                        s.setHWEP(Double.parseDouble(hwe));

                        s.setNrCalled(Integer.parseInt(nrCalled));
                    } catch (NumberFormatException e) {
//                    System.err.println("ERROR: cannot format number. ");
                    }
                    allSNPs.add(snp);
                    allSNPsInDs.add(s);
                    visitedSNPs.add(snp);
                }

                elems = tfIn.readLineElems(TextFile.tab);
            }
            tfIn.close();
            snps[dsN] = allSNPsInDs.toArray(new SNP[0]);

        }

        // hash the SNPs for easy lookup..
        snpIndex = new Integer[datasets.size()][allSNPs.size()]; // for each SNP in the meta analysis, this index points to a SNP in a dataset (in the jagged allSNPs array)
        snpToId = new HashMap<String, Integer>();
        int ctr = 0;
        for (String snp : allSNPs) {
            snpToId.put(snp, ctr);
            ctr++;
        }
        System.out.println(snpToId.size() + " snps loaded");
        for (int d = 0; d < datasets.size(); d++) {
            SNP[] snpsInDs = snps[d];
            for (int s = 0; s < snpsInDs.length; s++) {
                Integer id = snpToId.get(snps[d][s].getName());
                snpIndex[d][id] = s;
            }
        }

        // hash the eQTLs as well..
        eQTLIndex = new Integer[datasets.size()][allEQTL.size()];
        eQTLToId = new HashMap<String, Integer>();
        String[] allEQTLArr = allEQTL.toArray(new String[0]);
        for (int e = 0; e < allEQTLArr.length; e++) {
            String eQTL = allEQTLArr[e];
            eQTLToId.put(eQTL, e);
        }

        for (int d = 0; d < datasets.size(); d++) {
            List<String> eqtls = datasets.get(d).colObjects;
            for (int s = 0; s < eqtls.size(); s++) {
                String eQTL = eqtls.get(s);
                Integer id = null;
                if (isHT12v3[d]) {
                    id = eQTLToId.get(eQTL);
                } else {
                    String[] eqtlelems = eQTL.split("-");
                    eqtlelems[1] = ht12v4ToHT12v3.get(eqtlelems[1]);
                    eQTL = Strings.concat(eqtlelems, Strings.dash);
                    id = eQTLToId.get(eQTL);
                }
                if (id == null) {
                    System.err.println("ERROR: " + id + " is null for " + eQTL + " in dataset " + datasetFileDirs[d]);
                } else {
                    eQTLIndex[d][id] = s;
                }
            }
        }

        // hash the probes / covariates...
        allCovariates = new HashSet<String>();
        for (int d = 0; d < datasets.size(); d++) {
            DoubleMatrixDataset<String, String> ds = datasets.get(d);
            for (int r = 0; r < ds.nrRows; r++) {

                // each matrix has two 'special' covariates
                if (ds.rowObjects.get(r).equals("CellTypeInteractionZScore") || ds.rowObjects.get(r).equals("MainEffectZScore") || ds.rowObjects.get(r).equals("CellTypeZScore") || ds.rowObjects.get(r).equals("CellTypeSNPZScore")) {
                    allCovariates.add(ds.rowObjects.get(r));
                } else {
                    // parse the 'normal' covariates
                    if (!isHT12v3[d]) {
                        String cov = ht12v4ToHT12v3.get(ds.rowObjects.get(r));
                        if (cov == null) {
                            System.err.println("ERROR: HT12v4 probe not included in translation file? " + cov + "\t" + ds.rowObjects.get(r));
                            System.exit(-1);
                        }
                        if (cov.equals("-")) {
                            System.err.println("ERROR: covariate equals '-' " + ds.rowObjects.get(r));
                        } else {
                            allCovariates.add(cov);
                        }
                    } else {
                        String cov = ds.rowObjects.get(r);
                        if (cov.equals("-")) {
                            System.err.println("ERROR: covariate equals '-' " + ds.rowObjects.get(r));
                        } else {
                            allCovariates.add(cov);
                        }
                    }
                }
            }
        }

        allCovariateArr = allCovariates.toArray(new String[0]);
        covariateIndex = new Integer[datasets.size()][allCovariates.size()];
        covariateToId = new HashMap<String, Integer>();
        for (int c = 0; c < allCovariateArr.length; c++) {
            covariateToId.put(allCovariateArr[c], c);
        }

        for (int d = 0; d < datasets.size(); d++) {
            DoubleMatrixDataset<String, String> ds = datasets.get(d);
            for (int r = 0; r < ds.nrRows; r++) {
                String cov = null;

                if (ds.rowObjects.get(r).equals("CellTypeInteractionZScore") || ds.rowObjects.get(r).equals("MainEffectZScore") || ds.rowObjects.get(r).equals("CellTypeZScore") || ds.rowObjects.get(r).equals("CellTypeSNPZScore")) {
                    cov = ds.rowObjects.get(r);
                } else {
                    if (!isHT12v3[d]) {
                        cov = ht12v4ToHT12v3.get(ds.rowObjects.get(r));
                        if (cov == null) {
                            System.out.println("Covariate not found: " + ds.rowObjects.get(r));
                            System.exit(-1);
                        }
                    } else {
                        cov = ds.rowObjects.get(r);
                    }
                }

                if (cov != null) {
                    Integer id = covariateToId.get(cov);
                    if (id != null) {
                        covariateIndex[d][id] = r;
                    }
                }
            }
        }
    }

    public String[] getAllCovariateArr() {
        return allCovariateArr;
    }

    public Integer[][] getCovariateIndex() {
        return covariateIndex;
    }

    public Integer[][] geteQTLIndex() {
        return eQTLIndex;
    }

    public HashMap<String, Integer> geteQTLToId() {
        return eQTLToId;
    }

    public Integer[][] getSnpIndex() {
        return snpIndex;
    }

    public HashMap<String, Integer> getSnpToId() {
        return snpToId;
    }

    public HashSet<String> getAllSNPs() {
        return allSNPs;
    }

    public HashSet<String> getAllEQTL() {
        return allEQTL;
    }

    public ArrayList<DoubleMatrixDataset<String, String>> getDatasets() {
        return datasets;
    }

    public SNP[][] getSnps() {
        return snps;
    }

    public HashMap<String, Integer> getCovariateToId() {
        return covariateToId;
    }

    public HashSet<String> getAllCovariates() {
        return allCovariates;
    }
}
