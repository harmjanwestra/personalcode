package nl.harmjanwestra.playground.biogen.covariates;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import nl.harmjanwestra.playground.methylation.PCA;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.VIF;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class RewriteMDS {


    public static void main(String[] args) {

//		String[] linksfiles = new String[]{
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-MAYO-V2.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-MSBB-V2.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-ROSMAP-V2.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-Braineac.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-BrainGVEX.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-CMC.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-CMC_HBCC_set2.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-ENA.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-GTEx.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-GVEX.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-LIBD_1M.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-NABEC-H610.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-NABEC-WES.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-TargetALS.txt"
//		};

        String[] linksFilesEurCortex = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-MAYO-V2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-MSBB-V2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-ROSMAP-V2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-Bipseq_1M.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-Braineac.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-BrainGVEX-V2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-CMC.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-CMC_HBCC_set2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-ENA.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-GTEx.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-GVEX.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-Integrative-RemainingSamplesMergedIn.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-LIBD_1M.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-NABEC-H550.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-NABEC-H610.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-TargetALS.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-UCLA_ASD.txt"
        };

        String[] mdsfilesEurCortex = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MAYO-V2\\AMPAD-MAYO-V2-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MSBB-V2\\AMPAD-MSBB-V2-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-ROSMAP-V2\\AMPAD-ROSMAP-V2-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\Bipseq_1M\\Bipseq_1M-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\Braineac\\Braineac-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\BrainGVEX-V2\\BrainGVEXv2-mds.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC\\CMC-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC_HBCC_set2\\CMC_HBCC_set2-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\ENA\\ENA-postqc-mds-ibd.mds.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GVEX\\GVEX-mds-ibd.mds.gz",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\Integrative\\Integrative-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_1M\\LIBD_1M-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\NABEC-H550\\NABEC-H550-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\NABEC-H610\\NABEC-H610-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\TargetALS\\TargetALS-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\UCLA_ASD\\UCLA_ASD-mds-ibd.mds.gz"
        };

        String[] namesEurCortex = new String[]{
                "AMPAD-MAYO",
                "AMPAD-MSBB",
                "AMPAD-ROSMAP",
                "BipSeq-1M",
                "Braineac",
                "BrainGVEX-V2",
                "CMC",
                "CMC_HBCC_set2",
                "ENA",
                "GTEx",
                "GVEX",
//				"PsychEncode",
                "LIBD_1M",
                "NABEC-H550",
                "NABEC-H610",
                "TargetALS",
                "UCLA-ASD"
        };

        String[] linksFilesAFRCortex = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-AFR.txt-AMPAD-MSBB-V2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-AFR.txt-CMC.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-AFR.txt-CMC_HBCC_set1.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-AFR.txt-CMC_HBCC_set2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-AFR.txt-CMC_HBCC_set3.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-AFR.txt-GTEx.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-AFR.txt-LIBD_1M.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cortex.txt-dedup-gte.txt-AFR.txt-LIBD_h650.txt"
        };


        String[] mdsfilesAFRCortex = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MSBB-V2\\AMPAD-MSBB-V2-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC\\CMC-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC_HBCC_set1\\CMC_HBCC_set1-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC_HBCC_set2\\CMC_HBCC_set2-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC_HBCC_set3\\CMC_HBCC_set3-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_1M\\LIBD_1M-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_h650\\LIBD_h650-mds-ibd.mds.gz",
        };

        String[] namesAFRCortes = new String[]{
                "AMPAD-MSBB",
                "CMC",
                "CMC_HBCC_set1",
                "CMC_HBCC_set2",
                "CMC_HBCC_set3",
                "GTEx",
                "LIBD_1M",
                "LIBD_h650"
        };

        String outfileEurCortex = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-03-CovarsAndMDS\\2019-07-08-CovarsAndMDS-Cortex-EUR.txt";
        String outfileAFRCortex = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-03-CovarsAndMDS\\2019-07-01-CovarsAndMDS-Cortex-AFR.txt";


        String covars = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-11-Freeze2RNAQc\\2019-04-11-Freeze2.TMM.Covariates-Numeric-Top10Covariates.txt";

        RewriteMDS m = new RewriteMDS();

        String[] mdsfilesAll = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MAYO-V2\\AMPAD-MAYO-V2-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MSBB-V2\\AMPAD-MSBB-V2-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-ROSMAP-V2\\AMPAD-ROSMAP-V2-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\Bipseq_1M\\Bipseq_1M-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\Bipseq_h650\\Bipseq_h650-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\Braineac\\Braineac-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\BrainGVEX-V2\\BrainGVEXv2-mds.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC\\CMC-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC_HBCC_set1\\CMC_HBCC_set1-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC_HBCC_set2\\CMC_HBCC_set2-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC_HBCC_set3\\CMC_HBCC_set3-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\ENA\\ENA-postqc-mds-ibd.mds.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEX-WES\\GTEXWES-mds.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GVEX\\GVEX-mds-ibd.mds.gz",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\Integrative\\Integrative-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_1M\\LIBD_1M-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_h650\\LIBD_h650-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\NABEC-H550\\NABEC-H550-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\NABEC-H610\\NABEC-H610-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\NABEC-NeuroChip\\neurochip-mds.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\NABEC-WES\\nabec-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\TargetALS\\TargetALS-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\UCLA_ASD\\UCLA_ASD-mds-ibd.mds.gz"
        };

        String[] namesAll = new String[]{
                "AMPAD-MAYO-V2",
                "AMPAD-MSBB-V2",
                "AMPAD-ROSMAP-V2",
                "Bipseq_1M",
                "Bipseq_h650",
                "Braineac",
                "BrainGVEX-V2",
                "CMC",
                "CMC_HBCC_set1",
                "CMC_HBCC_set2",
                "CMC_HBCC_set3",
                "ENA",
                "GTEx",
                "GTEx-WES",
                "GVEX",
                "LIBD_1M",

                "LIBD_h650",
                "NABEC-H550",
                "NABEC-H610",
                "NABEC-NeuroChip",
                "NABEC-WES",
                "TargetALS",
                "UCLA_ASD"
        };

        String[] linkFilesAll = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\AMPAD-MAYO-V2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\AMPAD-MSBB-V2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\AMPAD-ROSMAP-V2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\Bipseq_1M.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\Bipseq_h650.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\Braineac.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\BrainGVEX-V2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\CMC.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\CMC_HBCC_set1.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\CMC_HBCC_set2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\CMC_HBCC_set3.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\ENA.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\GTEx.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\GTEx-WES.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\GVEX.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\LIBD_1M.txt",

                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\LIBD_h650.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\NABEC-H550.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\NABEC-H610.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\NABEC-NeuroChip.txt.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\NABEC-WES.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\TargetALS.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\allSplitPerDataset\\UCLA_ASD.txt"
        };

        String outfileall = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-03-CovarsAndMDS\\2019-07-01-CovarsAndMDS-All-All.txt";

        String[] linksFilesEURCerebellum = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cerebellum.txt-dedup-gte.txt-EUR.txt-AMPAD-MAYO-V2.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cerebellum.txt-dedup-gte.txt-EUR.txt-ENA.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cerebellum.txt-dedup-gte.txt-EUR.txt-GTEx.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cerebellum.txt-dedup-gte.txt-EUR.txt-TargetALS.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\splitpertissue\\Cerebellum.txt-dedup-gte.txt-EUR.txt-UCLA_ASD.txt"
        };
        String[] namesEURCerebellum = new String[]{
                "AMPAD-MAYO",
                "ENA",
                "GTEx",
                "TargetALS",
                "UCLA-ASD"
        };
        String[] mdsFilesEURCerebellum = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MAYO-V2\\AMPAD-MAYO-V2-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\ENA\\ENA-postqc-mds-ibd.mds.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\TargetALS\\TargetALS-mds-ibd.mds.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\UCLA_ASD\\UCLA_ASD-mds-ibd.mds.gz"
        };
        String outfileEurCerebellum = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-03-CovarsAndMDS\\2019-10-28-CovarsAndMDS-EUR-Cerebellum.txt";


        try {
//			m.run(linksFilesEurCortex, namesEurCortex, mdsfilesEurCortex, covars, outfileEurCortex);
//			m.run(linksFilesAFRCortex, namesAFRCortes, mdsfilesAFRCortex, covars, outfileAFRCortex);
//			m.run(linkFilesAll, namesAll, mdsfilesAll, covars, outfileall);
//			m.run(linksFilesEURCerebellum, namesEURCerebellum, mdsFilesEURCerebellum, covars, outfileEurCerebellum);


            // Freeze 2.1 - cortex
            String[] linkfiles = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-AFR.txt-AMPAD-MSBB-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-AFR.txt-CMC.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-AFR.txt-LIBD_1M.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-AFR.txt-LIBD_h650.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-MAYO-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-MSBB-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-ROSMAP-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-BrainGVEX-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-CMC.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-GTEx.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-GVEX.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-ENA.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-LIBD_1M.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-LIBD_h650.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-NABEC-H550.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-NABEC-H610.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-TargetALS.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-UCLA_ASD.txt"
            };

            String[] datasetnames = new String[]{
                    "AMPAD-MSBB-V2-AFR",
                    "CMC-AFR",
                    "LIBD_1M-AFR",
                    "LIBD_h650-AFR",
                    "AMPAD-MAYO-V2-EUR",
                    "AMPAD-MSBB-V2-EUR",
                    "AMPAD-ROSMAP-V2-EUR",
                    "BrainGVEX-V2-EUR",
                    "CMC-EUR",
                    "GTEx-EUR",
                    "GVEx",
                    "ENA-EUR",
                    "LIBD_1M-EUR",
                    "LIBD_h650-EUR",
                    "NABEC-H550-EUR",
                    "NABEC-H610-EUR",
                    "TargetALS-EUR",
                    "UCLA_ASD-EUR"
            };
            String[] mdsfiles = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MSBB-V2\\AMPAD-MSBB-V2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC\\CMC-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_1M\\LIBD_1M-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_h650\\LIBD_h650-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MAYO-V2\\AMPAD-MAYO-V2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MSBB-V2\\AMPAD-MSBB-V2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-ROSMAP-V2\\AMPAD-ROSMAP-V2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\BrainGVEX-V2\\BrainGVEXv2-mds.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC\\CMC-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GVEX\\GVEX-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-01-12-ENAGenotypes\\7-plinkMAF0001\\ENA-mds.mds",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_1M\\LIBD_1M-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_h650\\LIBD_h650-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\NABEC-H550\\NABEC-H550-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\NABEC-H610\\NABEC-H610-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\TargetALS\\TargetALS-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\UCLA_ASD\\UCLA_ASD-mds-ibd.mds.gz"
            };
            String covarfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-CovarsAndMDS\\2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated.txt.gz";
            String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-CovarsAndMDS\\2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cortex-withMDS.txt.gz";

            //   m.run(linkfiles, datasetnames, mdsfiles, covarfile, output);


            // Freeze 2.1 - cerebellum
            linkfiles = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cerebellum.txt-dedup-gte.txt-EUR.txt-AMPAD-MAYO-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cerebellum.txt-dedup-gte.txt-EUR.txt-GTEx.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cerebellum.txt-dedup-gte.txt-EUR.txt-TargetALS.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\cerebellum.txt-dedup-gte.txt-EUR.txt-UCLA_ASD.txt"
            };
            datasetnames = new String[]{
                    "AMPAD-MAYO-V2",
                    "GTEx",
                    "TargetALS",
                    "UCLA_ASD"
            };
            mdsfiles = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MAYO-V2\\AMPAD-MAYO-V2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\TargetALS\\TargetALS-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\UCLA_ASD\\UCLA_ASD-mds-ibd.mds.gz"
            };

            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-CovarsAndMDS\\2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cerebellum-withMDS.txt.gz";
            //   m.run(linkfiles, datasetnames, mdsfiles, covarfile, output);

            // Freeze 2.1 - basalganglia
            linkfiles = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\basalganglia.txt-dedup-gte.txt-EUR.txt-Braineac.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\basalganglia.txt-dedup-gte.txt-EUR.txt-GTEx.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\basalganglia.txt-dedup-gte.txt-EUR.txt-ENA.txt"
            };
            datasetnames = new String[]{
                    "Braineac",
                    "GTEx",
                    "ENA"
            };
            mdsfiles = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\Braineac\\Braineac-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-01-12-ENAGenotypes\\7-plinkMAF0001\\ENA-mds.mds"
            };

            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-CovarsAndMDS\\2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-basalganglia-withMDS.txt.gz";
            // m.run(linkfiles, datasetnames, mdsfiles, covarfile, output);


            // Freeze 2.1 - hippocampus
            linkfiles = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\hippocampus.txt-dedup-gte.txt-EUR.txt-Braineac.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output\\splitpertissuemorethan30samples\\hippocampus.txt-dedup-gte.txt-EUR.txt-GTEx.txt",
            };
            datasetnames = new String[]{
                    "Braineac",
                    "GTEx",
            };
            mdsfiles = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\Braineac\\Braineac-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
            };
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-CovarsAndMDS\\2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-hippocampus-withMDS.txt.gz";
//            m.run(linkfiles, datasetnames, mdsfiles, covarfile, output);

            ////// Freeze 2.1 population fix


            String[] tissues = new String[]{
                    "basalganglia",
                    "cerebellum",
                    "cortex",
                    "hippocampus",
                    "spinalcord"
            };

            String[] allLinkFiles = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\basalganglia.txt-dedup-gte.txt-EUR.txt-Braineac.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\basalganglia.txt-dedup-gte.txt-EUR.txt-GTEx.txt",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cerebellum.txt-dedup-gte.txt-EUR.txt-AMPAD-MAYO-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cerebellum.txt-dedup-gte.txt-EUR.txt-GTEx.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cerebellum.txt-dedup-gte.txt-EUR.txt-TargetALS.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cerebellum.txt-dedup-gte.txt-EUR.txt-UCLA_ASD.txt",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-AFR.txt-AMPAD-MSBB-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-AFR.txt-CMC_HBCC_set1.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-AFR.txt-CMC_HBCC_set2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-AFR.txt-CMC.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-AFR.txt-ENA.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-AFR.txt-LIBD_1M.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-AFR.txt-LIBD_h650.txt",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EAS.txt-ENA.txt",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-MAYO-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-MSBB-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-ROSMAP-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-BrainGVEX-V2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-CMC_HBCC_set2.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-CMC_HBCC_set3.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-CMC.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-ENA.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-GTEx.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-GVEX.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-LIBD_1M.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-LIBD_h650.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-NABEC-H550.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-NABEC-H610.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-TargetALS.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\cortex.txt-dedup-gte.txt-EUR.txt-UCLA_ASD.txt",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\hippocampus.txt-dedup-gte.txt-EUR.txt-Braineac.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\hippocampus.txt-dedup-gte.txt-EUR.txt-GTEx.txt",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\spinalcord.txt-dedup-gte.txt-EUR.txt-GTEx.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissuemorethan30samples\\spinalcord.txt-dedup-gte.txt-EUR.txt-TargetALS.txt"
            };

            String[] allMDSfiles = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\Braineac\\Braineac-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MAYO-V2\\AMPAD-MAYO-V2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\TargetALS\\TargetALS-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\UCLA_ASD\\UCLA_ASD-mds-ibd.mds.gz",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MSBB-V2\\AMPAD-MSBB-V2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC_HBCC_set1\\CMC_HBCC_set1-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC_HBCC_set2\\CMC_HBCC_set2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC\\CMC-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\ENA\\ENA-postqc-mds-ibd.mds.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_1M\\LIBD_1M-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_h650\\LIBD_h650-mds-ibd.mds.gz",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\ENA\\ENA-postqc-mds-ibd.mds.txt",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MAYO-V2\\AMPAD-MAYO-V2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-MSBB-V2\\AMPAD-MSBB-V2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\AMPAD-ROSMAP-V2\\AMPAD-ROSMAP-V2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\BrainGVEX-V2\\BrainGVEXv2-mds.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC_HBCC_set2\\CMC_HBCC_set2-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC_HBCC_set3\\CMC_HBCC_set3-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\CMC\\CMC-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\ENA\\ENA-postqc-mds-ibd.mds.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GVEX\\GVEX-mds-ibd.mds.gz",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_1M\\LIBD_1M-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\LIBD_h650\\LIBD_h650-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\NABEC-H550\\NABEC-H550-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\NABEC-H610\\NABEC-H610-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\TargetALS\\TargetALS-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\UCLA_ASD\\UCLA_ASD-mds-ibd.mds.gz",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\Braineac\\Braineac-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",

                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\qcfiles\\TargetALS\\TargetALS-mds-ibd.mds.gz",
            };

            if (allLinkFiles.length != allMDSfiles.length) {
                System.out.println(allLinkFiles.length + " link files, " + allMDSfiles.length + " MDS files");
                System.exit(-1);
            }

            covarfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-CovarsAndMDS\\2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated.txt.gz";
            for (String tissue : tissues) {
                ArrayList<String> selectedLinkFiles = new ArrayList<String>();
                ArrayList<String> selectedMDSFiles = new ArrayList<>();
                ArrayList<String> datasetNames = new ArrayList<>();

                for (int i = 0; i < allLinkFiles.length; i++) {

                    if (allLinkFiles[i].contains(tissue)) {
                        selectedLinkFiles.add(allLinkFiles[i]);
                        selectedMDSFiles.add(allMDSfiles[i]);
                        File f = new File(allLinkFiles[i]);
                        String filename = f.getName();

                        // cortex.txt-dedup-gte.txt-EUR.txt-GTEx.txt
                        filename = filename.replaceAll(tissue + ".txt-dedup-gte.txt-", "");
                        // EUR.txt-GTEx.txt
                        String[] filenameElems = filename.split("-");
                        String pop = filenameElems[0];
                        pop = pop.replaceAll(".txt", "");
                        filename = filename.replaceAll(pop + ".txt-", "");
                        String ds = filename.replaceAll(".txt", "");
                        ds = ds + "-" + pop;
                        datasetNames.add(ds);
                    }
                }


                output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-25-CovarsAndMDS\\2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-" + tissue + "-withMDS.txt.gz";
                m.run(selectedLinkFiles.toArray(new String[0]),
                        datasetNames.toArray(new String[0]),
                        selectedMDSFiles.toArray(new String[0]),
                        covarfile,
                        output);
            }


        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void run(String[] linkfiles, String[] datasetnames, String[] mdsfiles, String covarfile, String output) throws
            Exception {

        if (linkfiles.length != datasetnames.length) {
            System.out.println("Not equal lengths!");
            System.exit(-1);
        }

        DoubleMatrixDataset<String, String> covars = null;
        if (covarfile != null) {
            covars = DoubleMatrixDataset.loadDoubleData(covarfile);
        }

        int reference = 0;
        int nrlinks = 0;
        for (int i = 0; i < linkfiles.length; i++) {
            TextFile tf = new TextFile(linkfiles[i], TextFile.R);
            int size = tf.readAsArrayList().size();
            tf.close();
            if (size > nrlinks) {
                reference = i;
                nrlinks = size;

            }
        }
        System.out.println("Reference dataset is " + datasetnames[reference] + " with " + nrlinks + " links.");

        int[] datasets = new int[linkfiles.length - 1];
        TextFile tfout = new TextFile(output, TextFile.W);

        String header = "";
        if (covars != null) {
            String covarString = Strings.concat(covars.getColObjects(), Strings.tab);
            header = "Sample\t" + covarString + "\tMDS1\tMDS2\tMDS3\tMDS4";
        }


        for (int d = 0; d < datasetnames.length; d++) {
            if (d != reference) {
                header += "\t" + datasetnames[d];
            }
        }

        tfout.writeln(header);

        int total = 0;
        int totalcovars = 0;
        int dsctr = 0;
        for (int i = 0; i < linkfiles.length; i++) {
            datasets = new int[linkfiles.length - 1];
            HashMap<String, HashSet<String>> dnaToRNA = loadLinkFile(linkfiles[i]);
            if (i != reference) {
                datasets[dsctr] = 1;
                dsctr++;
            }

            int matched = 0;
            int hascovars = 0;
            // parse MDS file
            TextFile tf = new TextFile(mdsfiles[i], TextFile.R);
            tf.readLine(); // skip header

            System.out.println("Parsing: " + mdsfiles[i]);
            String line = tf.readLine();
            HashSet<String> foundActualDNAs = new HashSet<>();
            while (line != null) {

                while (line.startsWith(" ")) {
                    line = line.substring(1);
                }
                while (line.contains("  ")) {
                    line = line.replaceAll("  ", " ");
                }

                String[] elems = line.split(" ");
                String sample = elems[1];
                HashSet<String> idlist = dnaToRNA.get(sample);
                if (idlist == null) {
                    idlist = dnaToRNA.get("0_" + sample + "_" + sample);
                }


                boolean idhascovars = false;
                if (idlist == null && datasetnames[i].equals("ENA")) {
//					System.out.println("Could not find sample " + sample + " for dataset: " + datasetnames[i]);
                } else if (idlist != null) {

                    foundActualDNAs.add(sample);
                    for (String id : idlist) {

                        if (covars != null) {
                            Integer covarId = covars.getHashRows().get(id);
                            double[] covarvals = new double[covars.getColObjects().size()];
                            if (covarId != null) {
                                covars.getRow(covarId).toArray(covarvals);
                                idhascovars = true;
                            }
                            try {
                                String outln = id + "\t" +
                                        Strings.concat(covarvals, Strings.tab) + "\t" +
                                        Strings.concat(elems, Strings.tab, 3, 7) + "\t" +
                                        Strings.concat(datasets, Strings.tab);
                                tfout.writeln(outln);
                            } catch (ArrayIndexOutOfBoundsException e) {
                                System.out.println("Error: we don't appear to have all data for sample: " + id);
                                System.out.println("covars: " + covarvals.length);
                                System.out.println("elems: " + elems.length);
                                System.out.println("datasets: " + datasets.length);
                                System.exit(-1);
                            }
                        } else {
                            String outln = id + "\t" + Strings.concat(elems, Strings.tab, 3, 7) + "\t" + Strings.concat(datasets, Strings.tab);
                            tfout.writeln(outln);
                        }

                    }
                    matched++;
                    if (idhascovars) {
                        hascovars++;
                    }
                }
                line = tf.readLine();
            }
            totalcovars += hascovars;
            total += matched;

            tf.close();

            for (String s : dnaToRNA.keySet()) {
                if (!foundActualDNAs.contains(s)) {
                    System.out.println(datasetnames[i] + "\t" + s + "\tnot found.");
                }
            }

            System.out.println(datasetnames[i] + "\t" + matched + " matched out of " + dnaToRNA.size() + ", " + hascovars + " with covars");
        }

        tfout.close();

        System.out.println(total + "/" + totalcovars);

        System.out.println();
        System.out.println("Output is here: " + output);
        System.out.println();
        System.out.println("Now checking properties of new covars.");

        // remove non-variant columns
        DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(output);

//		boolean[] zerovariancecol = new boolean[ds.columns()];
//		int colsWithVariance = 0;
//		for (int c = 0; c < ds.columns(); c++) {
//			DoubleMatrix1D col = ds.getCol(c);
//			double var = Descriptives.variance(col.toArray());
//			if (var == 0) {
//				zerovariancecol[c] = true;
//			} else {
//				colsWithVariance++;
//			}
//		}
//		double[][] newMatrix = new double[colsWithVariance][ds.rows()];
//		ArrayList<String> newrows = new ArrayList<>();
//
//		int ctr = 0;
//		for (int c = 0; c < ds.columns(); c++) {
//			DoubleMatrix1D col = ds.getCol(c);
//			if (!zerovariancecol[c]) {
//				newMatrix[ctr] = col.toArray();
//				newrows.add(ds.getColObjects().get(c));
//				ctr++;
//			}
//		}

//		for (int c = 0; c < newMatrix.length; c++) {
//			double[] y = newMatrix[c];
//			double[][] x = new double[newMatrix.length - 1][];
//			int cctr = 0;
//			for (int c2 = 0; c2 < newMatrix.length; c2++) {
//				if (c2 != c) {
//					x[cctr] = newMatrix[c2];
//					cctr++;
//				}
//			}
//
//			omp.newSampleData(y, x);
//			double rsq = omp.calculateRSquared();
//			if (rsq >= 0.99){
//
//			}
//
//		}

        // check for collinearity
        VIF vif = new VIF();
        vif.vifCorrect(ds, 0.999);
        ds.save(output + "-VIFCorrected.txt.gz");

//

    }

    private HashMap<String, HashSet<String>> loadLinkFile(String linkfile) throws IOException {

        HashMap<String, HashSet<String>> dnaToRNA = new HashMap<String, HashSet<String>>();
        TextFile tf = new TextFile(linkfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);

        while (elems != null) {
            HashSet<String> list = dnaToRNA.get(elems[0]);
            if (list == null) {
                list = new HashSet<>();
            }
            list.add(elems[1]);
            dnaToRNA.put(elems[0], list);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        return dnaToRNA;

    }


}
