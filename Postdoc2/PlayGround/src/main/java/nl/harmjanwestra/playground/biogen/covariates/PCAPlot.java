package nl.harmjanwestra.playground.biogen.covariates;

import com.itextpdf.text.DocumentException;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.containers.Triple;
import umcg.genetica.graphics.*;
import umcg.genetica.graphics.panels.LegendPanel;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.graphics.themes.DefaultTheme;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class PCAPlot {

	public static void main(String[] args) {
//        String knownClasses = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p3v5a-sampleinfo.txt";
//        String superclasses = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p3v5a-superpopulations.txt";


//		String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\genotypepca\\plink.eigenvec.txt";
//		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\genotypepca\\AMPAD-";
//		String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\genotypepca\\plink.eigenvec";
//		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\genotypepca\\GTEx-";
//        String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\plink.eigenvec.txt";
//        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-";
//        String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\genotypepca\\plink.eigenvec.txt";
//        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\genotypepca\\CMC-";
//        String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\dnapca\\plink.eigenvec.txt";
//        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\dnapca\\Braineac-";


		PCAPlot p = new PCAPlot();

		try {
//            String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\hap550qc\\lng_4tissue_brain_qtl_fix-pca - Copy.eigenvec.txt";
//            String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\hap550-";
//
//            p.plot(pcafile, knownClasses, superclasses, output);
//            pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\hap660qc\\eQTL_v2_fix-pca.eigenvec";
//            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\hap610-";
//
//            p.plot(pcafile, knownClasses, superclasses, output);
//
//            String knownClasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-tissue.annot.txt";
//            String superclasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-tissue.annot-groups.txt";
//            String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\mergednorm\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.PCAOverSamplesEigenvectors.txt.gz";
////            String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\mergednorm\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors.txt.gz";
//            String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\MergedRNASeqPCAPlot-Centered-PC1v2-tissue.pdf";
//            p.plot(pcafile, knownClasses, superclasses, 1, 2, 7, output);
//            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\MergedRNASeqPCAPlot-Centered-PC3v4-tissue.pdf";
//            p.plot(pcafile, knownClasses, superclasses, 3, 4, 7, output);
//
//            knownClasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-dataset.annot.txt";
//            superclasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-dataset.annot-groups.txt";
//            pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\mergednorm\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.PCAOverSamplesEigenvectors.txt.gz";
////            String pcafile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\mergednorm\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors.txt.gz";
//            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\MergedRNASeqPCAPlot-Centered-PC1v2-dataset.pdf";
//            p.plot(pcafile, knownClasses, superclasses, 1, 2, 7, output);
//            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\MergedRNASeqPCAPlot-Centered-PC3v4-dataset.pdf";
//            p.plot(pcafile, knownClasses, superclasses, 3, 4, 7, output);


//            String inputxyfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcabeforecorrection\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.PCAOverSamplesEigenvectors.txt.gz";
//            String covariatefile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates-qualityscores-filter.txt";
//            String covariatefile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates-noncategorical.txt";
//            String outplot = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcabeforecorrection\\plotsbeforecorrection\\";
//            p.correlateCovariatesWithExp(inputxyfile, covariatefile, outplot, 0.5);

//			String inputxyfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcabeforecorrection\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.PCAOverSamplesEigenvectors.txt.gz";
////            String covariatefile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates-qualityscores-filter.txt";
//			String covariatefile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates-noncategorical.txt";
//			String outplot = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcabeforecorrection\\plotsbeforecorrection\\";
//			p.correlateCovariatesWithExp(inputxyfile, covariatefile, outplot, 0.5);


//
//            inputxyfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcadatasetlabel\\eigenvectorsAfterDatasetLabelCorrection.txt.gz";
//            outplot = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcadatasetlabel\\plots\\";
////            p.correlateCovariatesWithExp(inputxyfile, covariatefile, outplot, 0.5);
//
//            inputxyfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcacovariates\\EigenVectorsAfterCovarCorrection.txt.gz";
//            outplot = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcacovariates\\plots\\";
////            p.correlateCovariatesWithExp(inputxyfile, covariatefile, outplot, 0.5);
//
//            String knownClasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-dataset.annot.txt";
//            String superclasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-dataset.annot-groups.txt";
//
//            String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcacovariates\\plots\\MergedRNASeqPCAPlot-Centered-PC1v2-dataset.png";
//            p.plot(inputxyfile, knownClasses, superclasses, 1, 2, 7, output);
//            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcacovariates\\plots\\MergedRNASeqPCAPlot-Centered-PC3v4-dataset.png";
////            p.plot(inputxyfile, knownClasses, superclasses, 3, 4, 7, output);
//
//            knownClasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-tissue.annot.txt";
//            superclasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\rna-tissue.annot-groups.txt";
//            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcacovariates\\plots\\MergedRNASeqPCAPlot-Centered-PC1v2-tissue.png";
//            p.plot(inputxyfile, knownClasses, superclasses, 1, 2, 7, output);
//            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pcacovariates\\plots\\MergedRNASeqPCAPlot-Centered-PC3v4-tissue.png";
//            p.plot(inputxyfile, knownClasses, superclasses, 3, 4, 7, output);

			String knownClasses = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p3v5a-sampleinfo.txt";
			String superclasses = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p3v5a-superpopulations.txt";
			String[] eigenvecfiles = new String[]{
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\AMPAD-MAYO-V2\\AMPAD-MAYO-V2-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\AMPAD-MSBB-V2\\AMPAD-MSBB-V2-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\AMPAD-ROSMAP-V2\\AMPAD-ROSMAP-V2-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\Bipseq_1M\\Bipseq_1M-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\Bipseq_2pt5M\\Bipseq_2pt5M-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\Bipseq_5M\\Bipseq_5M-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\Bipseq_h650\\Bipseq_h650-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\Braineac\\Braineac-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\BrainGVEX\\BrainGVEX-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\CMC\\CMC-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\CMC_HBCC_set1\\CMC_HBCC_set1-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\CMC_HBCC_set2\\CMC_HBCC_set2-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\CMC_HBCC_set3\\CMC_HBCC_set3-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\GTEx\\GTEx-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\GVEX\\GVEX-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\GVEXUpdate\\GVEXUpdate-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\Integrative\\Integrative-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\LIBD_1M\\LIBD_1M-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\LIBD_5M\\LIBD_5M-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\LIBD_h650\\LIBD_h650-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\NABEC-H550\\NABEC-H550-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\NABEC-H610\\NABEC-H610-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\UCLA_ASD\\UCLA_ASD-pca.eigenvec.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\ENA\\ENA-vcf-pca.eigenvec",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\TargetALS\\TargetALS-pca.eigenvec",
					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\NABEC-WES\\NABEC-WES-pca.eigenvec.gz"
			};
			String[] outputs = new String[]{
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\AMPAD-MAYO-V2-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\AMPAD-MSBB-V2-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\AMPAD-ROSMAP-V2-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Bipseq_1M-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Bipseq_2pt5M-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Bipseq_5M-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Bipseq_h650-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Braineac-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\BrainGVEX-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\CMC-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\CMC_HBCC_set1-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\CMC_HBCC_set2-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\CMC_HBCC_set3-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\GTEx-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\GVEX-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\GVEXUpdate-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\Integrative-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\LIBD_1M-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\LIBD_5M-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\LIBD_h650-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\NABEC-H550-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\NABEC-H610-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\UCLA_ASD-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\ENA-vcf-",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\TargetALS-",
					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\populationassignment\\NABEC-WES-",

			};

//			for (int i = 0; i < eigenvecfiles.length; i++) {
//				String s = eigenvecfiles[i];
//				String output = outputs[i];
//				File path = new File(output);
//				String name = path.getName();
//				name = name.substring(0, name.length() - 1);
//				p.plot(s, knownClasses, superclasses, 2, 3, "1000 genomes", name, 7, output);
//			}
//			System.exit(0);
//
//			String inputxyfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\6-covcorrection\\pc1to4.txt";
//			knownClasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-brain.phenotypes.rna.txt";
//			superclasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-brain.phenotypes.rna-groups-CortexVsNonCortex.txt";
//			superclasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-brain.phenotypes.rna-groupsbiogen.txt";
//			String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\2019-04-16-Centered-CovCorrected-PC1v2-tissue-biogen.pdf";
//			p.plot(inputxyfile, knownClasses, superclasses, 1, 2, "Annotated Samples", "Unannotated samples", 7, output);
//			output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\2019-04-16-Centered-CovCorrected-PC3v4-tissue-biogen.pdf";
//			p.plot(inputxyfile, knownClasses, superclasses, 3, 4, "Annotated Samples", "Unannotated samples", 7, output);

//			String alltissues = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-brain.phenotypes.rna-groups2.txt";
//			output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\2019-04-16-Centered-CovCorrected-PC1v2-alltissue-test.png";

//			p.plot(inputxyfile, knownClasses, alltissues, 1, 2, "Annotated Samples", "Unannotated samples", 7, output);
//
//
//			knownClasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-brain.phenotypes.datasets.txt";
//			superclasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-brain.phenotypes.datasets-groups.txt";
//			output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\2019-04-13-Centered-CovCorrected-PC1v2-dataset.png";
//			p.plot(inputxyfile, knownClasses, superclasses, 1, 2, "Annotated Samples", "Unannotated samples", 7, output);
//			output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\2019-04-13-Centered-CovCorrected-PC3v4-dataset.png";
//			p.plot(inputxyfile, knownClasses, superclasses, 3, 4, "Annotated Samples", "Unannotated samples", 7, output);
//

			String inputxyfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\5-centerscalerun2\\5-pc1to4.txt";
//            String covariatefile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates-qualityscores-filter.txt";
			String covariatefile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-11-Freeze2.TMM.Covariates-Numeric-Top10Covariates.txt";
			String outplot = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\covcorrelation\\";
//			p.plotContinuousVariable(inputxyfile, covariatefile, outplot, 0.5);

//			p.plotColVsCol("D:\\Sync\\SyncThing\\Postdoc2\\2019-methylation\\pc1-10.txt", "D:\\Sync\\SyncThing\\Postdoc2\\2019-methylation\\pc1-10.png");


//			knownClasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-methylation\\outliers.txt";
//			superclasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-methylation\\GSEGroups.txt";
//			String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-methylation\\2019-05-09-PreCenterScale\\2019-05-09-32kPublicMethylationData-PreCenterScale-pc1-2-outliergroups.png";
//			inputxyfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-methylation\\2019-05-09-PreCenterScale\\pc1-10.txt";
//			p.plot(inputxyfile, knownClasses, superclasses, 1, 2, "outlier GSE", "other GSE", 7, output);
			knownClasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-methylation\\outliers.txt";
			superclasses = "D:\\Sync\\SyncThing\\Postdoc2\\2019-methylation\\GSEGroups.txt";
			String output = "D:\\TMP\\methpc\\pcs1-10-evs.png";
			String outputall = "D:\\TMP\\methpc\\pcs1-10-evs-allvsall.png";
			inputxyfile = "D:\\TMP\\methpc\\pcs1-10.txt";
			p.plotColVsCol(inputxyfile, outputall);
			p.plot(inputxyfile, knownClasses, superclasses, 1, 2, "outlier GSE", "other GSE", 7, output);

		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}


	private void plotColVsCol(String inputxyfile, String output) throws Exception {

		DefaultTheme def = new DefaultTheme();

		Color[] colors = new Color[]{

				Color.decode("#0099ff".toUpperCase()),
				Color.decode("#aaaaaa".toUpperCase()),
				Color.decode("#ffeb3b".toUpperCase()),
				Color.decode("#e91e63".toUpperCase()),
				Color.decode("#009688".toUpperCase()),
				Color.decode("#8bc34a".toUpperCase()),
				Color.decode("#ff5722".toUpperCase()),

				Color.decode("#3f51b5".toUpperCase()),
				Color.decode("#00bcd4".toUpperCase()),
				Color.decode("#9c27b0".toUpperCase()),
				Color.decode("#673ab7".toUpperCase()),

				Color.decode("#03a9f4".toUpperCase()),
				Color.decode("#607d8b".toUpperCase()),
				Color.decode("#4caf50".toUpperCase()),
				Color.decode("#cddc39".toUpperCase()),
				Color.decode("#ffc107".toUpperCase()),
				Color.decode("#ff9800".toUpperCase()),
				Color.decode("#795548".toUpperCase()),

		};

		int alpha = 250;
		for (int c = 0; c < colors.length; c++) {
			colors[c] = new Color(colors[c].getRed(), colors[c].getGreen(), colors[c].getBlue(), alpha);
		}

		def.setColors(colors);

		Color lg = def.getLightGrey();
		def.setLightgrey(new Color(lg.getRed(), lg.getGreen(), lg.getBlue(), alpha));

		DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(inputxyfile); // samples on rows


		Grid grid = new Grid(500, 500, ds.columns() + 1, ds.columns() + 1, 100, 100);


		for (int c = 0; c < ds.columns(); c++) {
			double[] colA = ds.getCol(c).toArray();
			for (int d = c + 1; d < ds.columns(); d++) {
				double[] colB = ds.getCol(d).toArray();
				ScatterplotPanel p = new ScatterplotPanel(1, 1);
				p.setData(colA, colB);
				p.setLabels(ds.getColObjects().get(c), ds.getColObjects().get(d));
				p.setAlpha(0.1f);
				p.setTheme(def);
				grid.addPanel(p, c, d);


			}
		}
		grid.draw(output, DefaultGraphics.Output.PNG);


	}

	private void plotContinuousVariable(String inputxyfile, String covariatefile, String outplot,
										double corthreshold) throws Exception {


		DefaultTheme def = new DefaultTheme();
		Color[] colors = new Color[]{
				Color.decode("#ff0000".toUpperCase()),
				Color.decode("#aaaaaa".toUpperCase()),
				Color.decode("#ffeb3b".toUpperCase()),
				Color.decode("#e91e63".toUpperCase()),
				Color.decode("#009688".toUpperCase()),
				Color.decode("#8bc34a".toUpperCase()),
				Color.decode("#ff5722".toUpperCase()),

				Color.decode("#3f51b5".toUpperCase()),
				Color.decode("#00bcd4".toUpperCase()),
				Color.decode("#9c27b0".toUpperCase()),
				Color.decode("#673ab7".toUpperCase()),

				Color.decode("#03a9f4".toUpperCase()),
				Color.decode("#607d8b".toUpperCase()),
				Color.decode("#4caf50".toUpperCase()),
				Color.decode("#cddc39".toUpperCase()),
				Color.decode("#ffc107".toUpperCase()),
				Color.decode("#ff9800".toUpperCase()),
				Color.decode("#795548".toUpperCase()),

		};

		int alpha = 120;
		for (int c = 0; c < colors.length; c++) {
			colors[c] = new Color(colors[c].getRed(), colors[c].getGreen(), colors[c].getBlue(), alpha);
		}

		def.setColors(colors);

		Color lg = def.getLightGrey();
		def.setLightgrey(new Color(lg.getRed(), lg.getGreen(), lg.getBlue(), alpha));

		DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(inputxyfile); // samples on rows
		DoubleMatrixDataset<String, String> dscovariate = DoubleMatrixDataset.loadDoubleData(covariatefile); // samples on rows


		int[] samplemap = new int[ds.rows()];
		for (int r = 0; r < ds.rows(); r++) {

			String rowname = ds.getRowObjects().get(r);

			Integer index = dscovariate.getHashRows().get(rowname);
			if (index == null) {
				index = -1;
			}
			samplemap[r] = index;
		}

		int nrpcs = 2;
		Grid grid = new Grid(500, 500, 3, 3, 100, 100);

		for (int c = 0; c < dscovariate.columns(); c++) {
			double[] covariate = dscovariate.getCol(c).toArray();

			String covariatename = dscovariate.getColObjects().get(c);
			double[] finalcovariate = new double[ds.rows()];
			double min = Double.MAX_VALUE;
			double max = -Double.MAX_VALUE;
			for (int i = 0; i < ds.rows(); i++) {
				int map = samplemap[i];
				if (map == -1) {
					finalcovariate[i] = Double.NaN;
				} else {
					double v = covariate[map];
					finalcovariate[i] = v;
					if (v < min) {
						min = v;
					}
					if (v > max) {
						max = v;
					}
				}
			}

			double sum = 0;
			int nrval = 0;
			double range = max - min;
			for (int d = 0; d < finalcovariate.length; d++) {
				double v = finalcovariate[d];
				if (!Double.isNaN(v)) {
					double perc = (v + min) / range;
					if (perc > 1) {
						perc = 1;
					}
					sum += perc;
					nrval++;
					finalcovariate[d] = perc;
				}
			}
			double avg = sum / nrval;
//            System.out.println("Error: " + v + "\t" + range + "\t" + min + "\t" + max + "\t" + perc);
//            System.out.print(Strings.concat(finalcovariate, Pattern.compile("\n")));
//            System.exit(-1);

			double maxcor = 0;
			DecimalFormat format = new DecimalFormat("#.###");
			for (int i = 0; i < nrpcs; i++) {
				double[] x = ds.getCol(i).toArray();
				double xcor = pruneAndCorrelate(x, finalcovariate);
				if (Math.abs(xcor) > maxcor) {
					maxcor = Math.abs(xcor);
				}
				for (int j = i + 1; j < nrpcs; j++) {
					double[] y = ds.getCol(j).toArray();
					double ycor = pruneAndCorrelate(y, finalcovariate);
					if (Math.abs(ycor) > maxcor) {
						maxcor = Math.abs(ycor);
					}
					ScatterplotPanel p = new ScatterplotPanel(1, 1);
					p.setTitle("x cor: " + format.format(xcor) + ", y cor: " + format.format(ycor));
					p.setLabels(ds.getColObjects().get(i), ds.getColObjects().get(j));
					p.setData(x, y);
					p.setColorValues(finalcovariate);
					grid.addPanel(p, i, j - 1);
					p.setTheme(def);
					p.setPlotElems(true, false);
				}
			}

			if (maxcor > corthreshold) {
				LegendPanel legend = new LegendPanel(1, 1);
				legend.setTheme(def);
				legend.setTitle(covariatename);
				legend.setVals(min, max);
				legend.setAvg(avg, nrval);
				grid.addPanel(legend, 1, 0);
				grid.draw(outplot + format.format(maxcor) + "-PCA-" + covariatename + ".pdf");
			}
		}


	}

	double pruneAndCorrelate(double[] x, double[] y) throws IOException, DocumentException {
		return pruneAndCorrelate(x, y, 1.1, null, null, null, true);
	}

	double pruneAndCorrelate(double[] x, double[] y, double corthreshold, String outputplot, String xlabel, String
			ylabel, boolean spearman) throws IOException, DocumentException {

		ArrayList<Double> xtmp = new ArrayList<>(x.length);
		ArrayList<Double> ytmp = new ArrayList<>(x.length);

		for (int i = 0; i < x.length; i++) {
			if (!Double.isNaN(x[i]) && !Double.isNaN(y[i])) {
				xtmp.add(x[i]);
				ytmp.add(y[i]);
			}
		}

		if (xtmp.size() > 3) {
			double[] xnew = Primitives.toPrimitiveArr(xtmp);
			double[] ynew = Primitives.toPrimitiveArr(ytmp);
			if (spearman) {
				SpearmansCorrelation c = new SpearmansCorrelation();
				double xcor = c.correlation(xnew, ynew);
				return xcor;
			} else {
				return Correlation.correlate(xnew, ynew);
			}
//			if (Math.abs(xcor) > corthreshold) {
//				System.out.println("Plot: " + xlabel + "\t" + ylabel + "\t" + xcor);
//				Grid g2 = new Grid(300, 300, 1, 1, 100, 100);
//				ScatterplotPanel sp = new ScatterplotPanel(1, 1);
//				sp.setData(xnew, ynew);
//				DecimalFormat f = new DecimalFormat("#.####");
//				sp.setTitle("r: " + f.format(xcor) + ", rsq: " + f.format(xcor * xcor) + ", n=" + xnew.length);
//				g2.addPanel(sp);
//				sp.setLabels(xlabel, ylabel);
//				sp.setPlotElems(true, false);
//				g2.draw(outputplot + xlabel + "-" + ylabel + ".png");
//			}
		} else {
			return 0;
		}


	}


	public void plot(String pcafile, String knownclasses, String superclasses, int col1, int col2, String
			referencename, String datasetname,
					 int startk, String output) throws IOException, DocumentException {

		// load super populations
		HashMap<String, String> superpop = new HashMap<>();
		TextFile tfc = new TextFile(superclasses, TextFile.R);
		tfc.readLine();
		String[] elemsf = tfc.readLineElems(TextFile.tab);
		while (elemsf != null) {
			superpop.put(elemsf[0], elemsf[2]);
			elemsf = tfc.readLineElems(TextFile.tab);
		}
		tfc.close();

		String nullpop = "NotDefined";


		// load known classes
		HashMap<String, HashSet<String>> samplesPerPop = new HashMap<>();
		HashMap<String, String> sampleToPop = new HashMap<>();
		TextFile tf = new TextFile(knownclasses, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			String sample = elems[0];
			String type = elems[1];
			type = superpop.get(type);
			HashSet<String> set = samplesPerPop.get(type);
			if (set == null) {
				set = new HashSet<>();
			}
			set.add(sample);
			samplesPerPop.put(type, set);
			sampleToPop.put(sample, type);

			elems = tf.readLineElems(Strings.whitespace);

		}
		tf.close();

		System.out.println(sampleToPop.size() + " with known population.");
		System.out.println(samplesPerPop.size() + " populations total.");


		ArrayList<String> populations = new ArrayList<>();
		ArrayList<String> populationstmp = new ArrayList<>();
		populationstmp.addAll(samplesPerPop.keySet());
		Collections.sort(populationstmp);
		populations.add(nullpop);
		populations.addAll(populationstmp);


		HashMap<String, Integer> populationIndex = new HashMap<>();
		ArrayList<ArrayList<Double>> xvals = new ArrayList<>();
		ArrayList<ArrayList<Double>> yvals = new ArrayList<>();
		ArrayList<ArrayList<Double>> xvalsunlabeled = new ArrayList<>();
		ArrayList<ArrayList<Double>> yvalsunlabeled = new ArrayList<>();
		for (int i = 0; i < populations.size(); i++) {
			String pop = populations.get(i);
			populationIndex.put(pop, i);
			xvals.add(new ArrayList<>());
			xvalsunlabeled.add(new ArrayList<>());
			yvals.add(new ArrayList<>());
			yvalsunlabeled.add(new ArrayList<>());
		}


		// parse pca file
		ArrayList<Triple<String, Double, Double>> data = new ArrayList<Triple<String, Double, Double>>();
		TextFile tf2 = new TextFile(pcafile, TextFile.R);

		String ln = tf2.readLine();
		double minX = Double.MAX_VALUE;
		double maxX = -Double.MAX_VALUE;
		double minY = Double.MAX_VALUE;
		double maxY = -Double.MAX_VALUE;
		int ctr = 0;
		while (ln != null) {
			String pcadata = ln;
			while (pcadata.contains("  ")) {
				ln.replaceAll("  ", " ");
			}
			String[] pcaelems = Strings.whitespace.split(pcadata);
			String sample = pcaelems[0];
			if (pcafile.endsWith("eigenvec") || pcafile.endsWith("eigenvec.gz")) {
				sample = pcaelems[1];
			}


			try {
				Double pca1 = Double.parseDouble(pcaelems[col1]);

				Double pca2 = Double.parseDouble(pcaelems[col2]);
				if (!Double.isNaN(pca1)) {
					if (pca1 > maxX) {
						maxX = pca1;
					}
					if (pca1 < minX) {
						minX = pca1;
					}
				}
				if (!Double.isNaN(pca2)) {
					if (pca2 > maxY) {
						maxY = pca2;
					}
					if (pca2 < minY) {
						minY = pca2;
					}
				}

				data.add(new Triple<>(sample, pca1, pca2));

				String population = sampleToPop.get(sample);
				if (population == null) {
					population = nullpop;
				}
				sampleToPop.put(sample, population);
				Integer popindex = populationIndex.get(population);
				if (population.equals(nullpop)) {
					xvalsunlabeled.get(popindex).add(pca1);
					yvalsunlabeled.get(popindex).add(pca2);
				} else {
					xvals.get(popindex).add(pca1);
					yvals.get(popindex).add(pca2);
				}
			} catch (NumberFormatException e) {
				System.out.println("Error parsing line: " + e.getMessage());
				System.out.println("Probably a header I didnt expect.");
			}
			ln = tf2.readLine();
			ctr++;
			if (ctr % 1000 == 0) {
				System.out.println(ctr + " lines parsed");
			}
		}
		tf2.close();


		Range r = new Range(minX, minY, maxX, maxY);
		DefaultTheme def = new DefaultTheme();
		Color[] colors = new Color[]{
				new Color(0, 0, 0),
				new Color(145, 145, 145),

				Color.decode("#f44336".toUpperCase()),
				Color.decode("#2196f3".toUpperCase()),
				Color.decode("#ffeb3b".toUpperCase()),
				Color.decode("#e91e63".toUpperCase()),
				Color.decode("#009688".toUpperCase()),
				Color.decode("#8bc34a".toUpperCase()),
				Color.decode("#ff5722".toUpperCase()),

				Color.decode("#3f51b5".toUpperCase()),
				Color.decode("#00bcd4".toUpperCase()),
				Color.decode("#9c27b0".toUpperCase()),
				Color.decode("#673ab7".toUpperCase()),

				Color.decode("#03a9f4".toUpperCase()),
				Color.decode("#607d8b".toUpperCase()),
				Color.decode("#4caf50".toUpperCase()),
				Color.decode("#cddc39".toUpperCase()),
				Color.decode("#ffc107".toUpperCase()),
				Color.decode("#ff9800".toUpperCase()),
				Color.decode("#795548".toUpperCase()),

		};
		def.setColors(colors);

		// scatterplot
		ScatterplotPanel p1 = preparePanel(r, populations, xvals, yvals, referencename, "PC" + col1, "PC" + col2, def);
		ScatterplotPanel p2 = preparePanel(r, populations, xvalsunlabeled, yvalsunlabeled, "Dataset", "PC" + col1, "PC" + col2, def);

		Grid grid = new Grid(500, 500, 2, 3, 100, 100);
		grid.addPanel(p1);
//		grid.addPanel(p2);

		// knn population assignment
		HashMap<String, String> sampleToPopTmp = assignPopulation(data, sampleToPop, nullpop, startk);

		xvalsunlabeled = new ArrayList<>();
		yvalsunlabeled = new ArrayList<>();
		for (int i = 0; i < populations.size(); i++) {
			xvalsunlabeled.add(new ArrayList<>());
			yvalsunlabeled.add(new ArrayList<>());
		}

		TextFile assignmentout = new TextFile(output + "SampleAssignment.txt", TextFile.W);
		TextFile assignmentout2 = new TextFile(output + "SampleAssignmentKnown.txt", TextFile.W);

		for (Triple<String, Double, Double> sample : data) {
			String prevpop = sampleToPop.get(sample.getLeft());
			if (prevpop.equals(nullpop)) {
				String pop = sampleToPopTmp.get(sample.getLeft());
				Integer popindex = populationIndex.get(pop);
				if (popindex == null) {
					System.out.println("Could not find population : " + pop + " for sample: " + sample.getLeft());
				} else {
					xvalsunlabeled.get(popindex).add(sample.getMiddle());
					yvalsunlabeled.get(popindex).add(sample.getRight());
					assignmentout.writeln(sample.getLeft() + "\t" + pop);
				}
			} else {
				String pop = sampleToPopTmp.get(sample.getLeft());
				Integer popindex = populationIndex.get(pop);
				if (popindex == null) {
					System.out.println("Could not find population : " + pop + " for sample: " + sample.getLeft());
				} else {
					assignmentout2.writeln(sample.getLeft() + "\t" + pop);
				}
			}
		}
//        for (Triple<String, Double, Double> sample : data2) {
//            String prevpop = sampleToPop.get(sample.getLeft());
//            if (prevpop.equals(nullpop)) {
//                String pop = sampleToPopTmp.get(sample.getLeft());
//                Integer popindex = populationIndex.get(pop);
//                if (popindex == null) {
//                    System.out.println("Could not find population : " + pop + " for sample: " + sample.getLeft());
//                } else {
//                    assignmentout.writeln(sample.getLeft() + "\t" + pop);
//                }
//            }
//        }
		assignmentout2.close();
		assignmentout.close();

		ScatterplotPanel p3 = preparePanel(r, populations, xvalsunlabeled, yvalsunlabeled, datasetname, "PC" + col1, "PC" + col2, def);
		grid.addPanel(p3);


		grid.draw(output + "SampleAssignment.pdf");
	}

	private ScatterplotPanel preparePanel(Range datarange,
										  ArrayList<String> populations,
										  ArrayList<ArrayList<Double>> x,
										  ArrayList<ArrayList<Double>> y,
										  String title,
										  String xlabel,
										  String ylabel,
										  DefaultTheme def) {
		double[][] xprim = new double[populations.size()][];
		double[][] yprim = new double[populations.size()][];
		for (int i = 0; i < x.size(); i++) {
			xprim[i] = Primitives.toPrimitiveArr(x.get(i));
			yprim[i] = Primitives.toPrimitiveArr(y.get(i));
		}

		ScatterplotPanel punlab = new ScatterplotPanel(1, 1);
		punlab.setDataRange(datarange);
		punlab.setData(xprim, yprim);
		punlab.setTitle(title);
		punlab.setDatasetLabels(populations.toArray(new String[0]));
		punlab.setLabels(xlabel, ylabel);
		punlab.setPlotElems(true, true);
		punlab.setTheme(def);
		punlab.setAlpha(0.8f);

		return punlab;
	}

	private HashMap<String, String> assignPopulation
			(ArrayList<Triple<String, Double, Double>> data, HashMap<String, String> sampleToPop, String nullpop,
			 int startk) {
		int samplesWithoutPop = 0;
		HashMap<String, String> sampleToPopTmp = new HashMap<>();

		for (Triple<String, Double, Double> sample : data) {
			String pop = sampleToPop.get(sample.getLeft());
			if (pop.equals(nullpop)) {
				int k = startk;
				boolean resolved = false;
				while (!resolved) {
					String[] neighbors = new String[k];
					double[] distances = new double[k];
					double maxdist = -Double.MAX_VALUE;
					int n = 0;


					// this sample needs assignment
					// get startk nearest neighbors that do have an assignment
					for (Triple<String, Double, Double> sample2 : data) {
						String pop2 = sampleToPop.get(sample2.getLeft());
						if (!pop2.equals(nullpop)) {
							double distance = Math.sqrt(sq(sample.getMiddle() - sample2.getMiddle()) + sq(sample.getRight() - sample2.getRight()));

							if (n < k) {
								neighbors[n] = sample2.getLeft();
								distances[n] = distance;
								if (distance > maxdist) {
									maxdist = distance;
								}
								n++;
							} else if (distance < maxdist) {
								boolean update = false;
								double newmax = -Double.MAX_VALUE;

								for (int i = 0; i < neighbors.length; i++) {
									if (!update && (distances[i] > distance)) {
										distances[i] = distance;
										neighbors[i] = sample2.getLeft();
										update = true;
									}

									if (distances[i] > newmax) {
										newmax = distances[i];
									}
								}
								maxdist = newmax;
							}
						}
					}

					HashMap<String, Integer> ctr = new HashMap<>();
					HashMap<String, Double> sumd = new HashMap<>();

					for (int i = 0; i < neighbors.length; i++) {
						String pop2 = sampleToPop.get(neighbors[i]);
						Integer ct = ctr.get(pop2);
						Double sum = sumd.get(pop2);
						if (ct == null) {
							ct = 1;
							sum = distances[i];
						} else {
							ct++;
							sum += distances[i];
						}
						ctr.put(pop2, ct);
						sumd.put(pop2, sum);
					}


					int maxct = 0;
					String maxpop = null;

					for (String key : ctr.keySet()) {
						Integer ct = ctr.get(key);
						if (ct > maxct) {
							maxpop = key;
							maxct = ct;
						} else if (ct == maxct) {
							// break ties using absolute sum of distance
							double sum = sumd.get(key);
							if (sum < sumd.get(maxpop)) {
								maxpop = key;
								maxct = ct;
							}
						}
					}

					// assign
					sampleToPopTmp.put(sample.getLeft(), maxpop);
					resolved = true;
				}
			} else {

				sampleToPopTmp.put(sample.getLeft(), pop);
			}
		}

		return sampleToPopTmp;
	}

	private double sq(double v) {
		return v * v;
	}


}
