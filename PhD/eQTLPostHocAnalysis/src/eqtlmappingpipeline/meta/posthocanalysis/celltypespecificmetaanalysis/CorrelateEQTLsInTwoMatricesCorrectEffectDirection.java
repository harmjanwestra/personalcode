/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class CorrelateEQTLsInTwoMatricesCorrectEffectDirection {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here
            String dir1 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/EGCUT/";
            String dir2 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/Groningen/";
            String dir3 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/Rotterdam/";
            String out = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/EGCUTVsGroningenCorrelationOfInteractionTerms.txt";
            String pbt = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";

            CorrelateEQTLsInTwoMatricesCorrectEffectDirection c = new CorrelateEQTLsInTwoMatricesCorrectEffectDirection();
            boolean[] isHt12v4 = new boolean[]{false, false, true};
            c.run(dir1, dir2, dir3, pbt, isHt12v4, out);
        } catch (IOException ex) {
            Logger.getLogger(CorrelateEQTLsInTwoMatricesCorrectEffectDirection.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
    HashSet<String> idsToAvoid = new HashSet<String>();

    public void run(String dir1, String dir2, String dir3, String probetranslationfile, boolean[] isht12v4, String output) throws IOException {
        DoubleMatrixDataset<String, String> matrix1 = new DoubleMatrixDataset<String, String>(dir1 + "CellTypeSpecificityMatrix.binary");
        DoubleMatrixDataset<String, String> matrix2 = new DoubleMatrixDataset<String, String>(dir2 + "CellTypeSpecificityMatrix.binary");
        DoubleMatrixDataset<String, String> matrix3 = new DoubleMatrixDataset<String, String>(dir3 + "CellTypeSpecificityMatrix.binary");

        String ht12v3String = "HumanHT-12_V3_0_R2_11283641_A.txt";
        String ht12v4String = "HumanHT-12_V4_0_R1_15002873_B.txt";

        String hugoStr = "HUGO";

        //       PROBE LOOKUP TABLES..
        ProbeTranslation pbt = new ProbeTranslation();
        HashMap<String, String> ht12v3ToHt12v4 = new HashMap<String, String>();
        HashMap<String, String> ht12v3ToHUGO = new HashMap<String, String>();
        ht12v3ToHt12v4 = pbt.getProbeTranslation(probetranslationfile, ht12v3String, ht12v4String);
        ht12v3ToHUGO = pbt.getProbeTranslation(probetranslationfile, ht12v3String, hugoStr);

        idsToAvoid.add("CellTypeZScore");
        idsToAvoid.add("MainEffectZScore");
        idsToAvoid.add("CellTypeInteractionZScore");
        idsToAvoid.add("CellTypeSNPZScore");


        SNP[] snpsInDir1 = loadSNPs(dir1);
        SNP[] snpsInDir2 = loadSNPs(dir2);
        SNP[] snpsInDir3 = loadSNPs(dir3);

        HashMap<String, SNP> snpToSNPId1 = new HashMap<String, SNP>();
        for (int i = 0; i < snpsInDir1.length; i++) {
            snpToSNPId1.put(snpsInDir1[i].getName(), snpsInDir1[i]);
        }

        HashMap<String, SNP> snpToSNPId2 = new HashMap<String, SNP>();
        for (int i = 0; i < snpsInDir2.length; i++) {
            snpToSNPId2.put(snpsInDir2[i].getName(), snpsInDir2[i]);
        }

        HashMap<String, SNP> snpToSNPId3 = new HashMap<String, SNP>();
        for (int i = 0; i < snpsInDir3.length; i++) {
            snpToSNPId3.put(snpsInDir3[i].getName(), snpsInDir3[i]);
        }



        matrix1.transposeDataset();
        matrix2.transposeDataset();
        matrix3.transposeDataset();
        
        Integer col1 = matrix1.hashCols.get("MainEffectZScore");
        Integer col2 = matrix2.hashCols.get("MainEffectZScore");
        Integer col3 = matrix3.hashCols.get("MainEffectZScore");
        

        TextFile tf = new TextFile(output, TextFile.W);
        for (int row = 0; row < matrix1.nrRows; row++) {
            String eqtl = matrix1.rowObjects.get(row);
            Integer rowIdInOtherMatrix2 = matrix2.hashRows.get(eqtl);
            Integer rowIdInOtherMatrix3 = null;
            String snp = eqtl.split("-")[0];
            String probe = eqtl.split("-")[1];

            String ht12v4 = ht12v3ToHt12v4.get(probe);
            if (ht12v4 != null) {
                rowIdInOtherMatrix3 = matrix3.hashRows.get(snp + "-" + ht12v4);
            }

            if (rowIdInOtherMatrix3 == null) {
                System.err.println(eqtl+"\t"+snp + "-" + ht12v4+" not found in ht12v4 dataset");
            }


            SNP snp1 = snpToSNPId1.get(snp);
            SNP snp2 = snpToSNPId2.get(snp);
            SNP snp3 = snpToSNPId3.get(snp);
            if (snp1 != null && snp2 != null && snp3 != null && rowIdInOtherMatrix2 != null && rowIdInOtherMatrix3 != null) {


                Boolean flip2 = BaseAnnot.flipalleles(BaseAnnot.getAllelesDescription(snp1.getAlleles()), BaseAnnot.toString(snp1.getMinorAllele()), BaseAnnot.getAllelesDescription(snp2.getAlleles()), BaseAnnot.toString(snp2.getMinorAllele()));
                Boolean flip3 = BaseAnnot.flipalleles(BaseAnnot.getAllelesDescription(snp1.getAlleles()), BaseAnnot.toString(snp1.getMinorAllele()), BaseAnnot.getAllelesDescription(snp3.getAlleles()), BaseAnnot.toString(snp3.getMinorAllele()));

                ArrayList<Double> x = new ArrayList<Double>();
                ArrayList<Double> y = new ArrayList<Double>();
                ArrayList<Double> z = new ArrayList<Double>();


                for (int col = 0; col < matrix1.nrCols; col++) {
                    String colName = matrix1.colObjects.get(col);
                    if (!idsToAvoid.contains(colName)) {
                        Integer idInOtherMatrix2 = matrix2.hashCols.get(colName);

                        ht12v4 = ht12v3ToHt12v4.get(colName);

                        Integer idInOtherMatrix3 = matrix3.hashCols.get(ht12v4);
                        if (idInOtherMatrix2 != null && idInOtherMatrix3 != null) {
                            x.add(matrix1.rawData[row][col]);
                            if (flip2 != null && flip2) {
                                y.add(-matrix2.rawData[rowIdInOtherMatrix2][idInOtherMatrix2]);
                            } else {
                                y.add(matrix2.rawData[rowIdInOtherMatrix2][idInOtherMatrix2]);
                            }

                            if (flip3 != null && flip3) {
                                z.add(-matrix3.rawData[rowIdInOtherMatrix3][idInOtherMatrix3]);
                            } else {
                                z.add(matrix3.rawData[rowIdInOtherMatrix3][idInOtherMatrix3]);
                            }
                        }
                    }


                }
                double[] xarr = Primitives.toPrimitiveArr(x.toArray(new Double[0]));
                double[] yarr = Primitives.toPrimitiveArr(y.toArray(new Double[0]));
                double[] zarr = Primitives.toPrimitiveArr(z.toArray(new Double[0]));
                double corr2 = JSci.maths.ArrayMath.correlation(xarr, yarr);
                double corr3 = JSci.maths.ArrayMath.correlation(xarr, zarr);
                
                double m2 = matrix2.rawData[rowIdInOtherMatrix2][col2];
                if(flip2){
                    m2*=-1;
                }
                
                double m3 = matrix3.rawData[rowIdInOtherMatrix3][col3];
                if(flip3){
                    m3*=-1;
                }
                
                String outputStr = eqtl + "\t" + flip2 + "\t" + xarr.length + "\t" + corr2 + "\t" + corr3 + "\t" + matrix1.rawData[row][col1] + "\t" + m2 + "\t" + m3;
                System.out.println(outputStr);
                tf.writeln(outputStr);
            }

        }

        tf.close();


    }

    private SNP[] loadSNPs(String dir) throws IOException {
        TextFile tfIn = new TextFile(dir + "SNPSummaryStatistics.txt", TextFile.R);
        tfIn.readLine(); // skip the header..
        ArrayList<SNP> allSNPsInDs = new ArrayList<SNP>();
        String[] elems = tfIn.readLineElems(TextFile.tab);
        while (elems != null) {
            // SNP	Chr	ChrPos	Alleles	MinorAllele	MAF	CallRate	HWE	GenotypesCalled
            String snp = elems[0];
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

            s.setMAF(Double.parseDouble(maf));
            s.setCR(Double.parseDouble(cr));
            s.setHWEP(Double.parseDouble(hwe));

            s.setNrCalled(Integer.parseInt(nrCalled));
            allSNPsInDs.add(s);
            elems = tfIn.readLineElems(TextFile.tab);
        }
        tfIn.close();
        return allSNPsInDs.toArray(new SNP[0]);
    }
}
