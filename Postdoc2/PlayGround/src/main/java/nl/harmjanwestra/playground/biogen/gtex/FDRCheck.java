package nl.harmjanwestra.playground.biogen.gtex;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.HistogramPanel;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class FDRCheck {
	
	public static void main(String[] args) {
		String eqtlgen = "D:\\eQTLTest\\2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved.txt.gz";
		String[] files = new String[]{
				"D:\\eQTLTest\\psychencode-DER-08a_hg19_eQTL.significant-withalleles.txt",
				"D:\\eQTLTest\\DER-08b_hg19_eQTL.bonferroni-withalleles.txt",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Adipose_Subcutaneous.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Adipose_Visceral_Omentum.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Adrenal_Gland.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Artery_Aorta.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Artery_Coronary.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Artery_Tibial.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Amygdala.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Anterior_cingulate_cortex_BA24.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Caudate_basal_ganglia.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Cerebellar_Hemisphere.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Cerebellum.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Cortex.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Frontal_Cortex_BA9.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Hippocampus.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Hypothalamus.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Nucleus_accumbens_basal_ganglia.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Putamen_basal_ganglia.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Spinal_cord_cervical_c-1.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Substantia_nigra.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Breast_Mammary_Tissue.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Cells_EBV-transformed_lymphocytes.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Cells_Transformed_fibroblasts.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Colon_Sigmoid.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Colon_Transverse.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Esophagus_Gastroesophageal_Junction.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Esophagus_Mucosa.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Esophagus_Muscularis.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Heart_Atrial_Appendage.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Heart_Left_Ventricle.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Liver.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Lung.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Minor_Salivary_Gland.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Muscle_Skeletal.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Nerve_Tibial.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Ovary.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Pancreas.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Pituitary.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Prostate.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Skin_Not_Sun_Exposed_Suprapubic.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Skin_Sun_Exposed_Lower_leg.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Small_Intestine_Terminal_Ileum.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Spleen.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Stomach.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Testis.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Thyroid.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Uterus.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Vagina.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Whole_Blood.v7.egenes.txt.gz",
			
		};
		
		String[] files2 = new String[]{
				"D:\\eQTLTest\\psychencode-DER-08a_hg19_eQTL.significant-withalleles.txt",
				"D:\\eQTLTest\\psychencode-DER-08b_hg19_eQTL.bonferroni-withalleles.txt",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Adipose_Subcutaneous.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Adipose_Visceral_Omentum.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Adrenal_Gland.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Artery_Aorta.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Artery_Coronary.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Artery_Tibial.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Amygdala.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Anterior_cingulate_cortex_BA24.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Caudate_basal_ganglia.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Cerebellar_Hemisphere.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Cerebellum.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Cortex.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Frontal_Cortex_BA9.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Hippocampus.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Hypothalamus.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Nucleus_accumbens_basal_ganglia.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Putamen_basal_ganglia.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Spinal_cord_cervical_c-1.v7.egenes.txt.gz",
				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Substantia_nigra.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Breast_Mammary_Tissue.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Cells_EBV-transformed_lymphocytes.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Cells_Transformed_fibroblasts.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Colon_Sigmoid.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Colon_Transverse.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Esophagus_Gastroesophageal_Junction.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Esophagus_Mucosa.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Esophagus_Muscularis.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Heart_Atrial_Appendage.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Heart_Left_Ventricle.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Liver.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Lung.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Minor_Salivary_Gland.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Muscle_Skeletal.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Nerve_Tibial.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Ovary.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Pancreas.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Pituitary.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Prostate.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Skin_Not_Sun_Exposed_Suprapubic.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Skin_Sun_Exposed_Lower_leg.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Small_Intestine_Terminal_Ileum.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Spleen.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Stomach.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Testis.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Thyroid.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Uterus.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Vagina.v7.egenes.txt.gz",
//				"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Whole_Blood.v7.egenes.txt.gz"
			
		};
		
		FDRCheck fdr = new FDRCheck();
		String out = "D:\\eQTLTest\\significanteqtls-replineqtlgen-GTEXv7-brainForPsychencode.txt";
		try {
//			fdr.run(files, out);
			files2 = new String[]{
//					"D:\\biogen\\metabrainfx\\cis\\eQTLProbesFDR0.05-ProbeLevel.txt.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLProbesFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz",
					"D:\\eQTLTest\\psychencode-DER-08a_hg19_eQTL.significant-withalleles.txt",
					"D:\\eQTLTest\\psychencode-DER-08b_hg19_eQTL.bonferroni-withalleles.txt",
//					"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Amygdala.v7.egenes.txt.gz",
//					"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Anterior_cingulate_cortex_BA24.v7.egenes.txt.gz",
//					"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Caudate_basal_ganglia.v7.egenes.txt.gz",
			};
			fdr.runCheckGTExInEQTLGen(files2, eqtlgen, out);
//			fdr.filtergenes("D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLProbesFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz");
//			String annotation = "";
			String eqtlgenannot = "D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
			String metabrainannot = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\annotation\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz";
			files2 = new String[]{
//					"D:\\biogen\\metabrainfx\\cis\\eQTLProbesFDR0.05-ProbeLevel.txt.gz",
					"D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLProbesFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz",
					"D:\\eQTLTest\\psychencode-DER-08a_hg19_eQTL.significant-withalleles.txt",
					"D:\\eQTLTest\\psychencode-DER-08b_hg19_eQTL.bonferroni-withalleles.txt",
//					"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Amygdala.v7.egenes.txt.gz",
//					"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Anterior_cingulate_cortex_BA24.v7.egenes.txt.gz",
//					"D:\\eQTLTest\\GTEx_Analysis_v7_eQTL\\Brain_Caudate_basal_ganglia.v7.egenes.txt.gz",
			};
//			fdr.tssdistance(files2, eqtlgenannot, metabrainannot, "D:\\eQTLTest\\significanteqtls-replineqtlgen-GTEXv7-tssdist-brainForPsychencode.pdf");

//			files2 = new String[]{
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLProbesFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz",
//					"D:\\eQTLTest\\psychencode-DER-08a_hg19_eQTL.significant-withalleles.txt",
//			};
//			fdr.tssdistancePerFDRBin(files2, eqtlgenannot, metabrainannot, "D:\\eQTLTest\\significanteqtls-replineqtlgen-GTEXv7-tssdist-brainForPsychencode-perFDRBin.pdf");
//			fdr.compareTSSDistance(files, annotation, out);
		} catch (Exception e) {
			e.printStackTrace();
		}
//		catch (DocumentException e) {
//			e.printStackTrace();
//		}
		
		
	}
	
	public void filtergenes(String in, String out) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfo = new TextFile(out, TextFile.W);
		tfo.writeln(tf.readLine());
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> genes = new HashSet<>();
		while (elems != null) {
			String gene = elems[4];
			if (!genes.contains(gene)) {
				tfo.writeln(Strings.concat(elems, Strings.tab));
				genes.add(gene);
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfo.close();
	}
	
	public void run(String[] files, String out) throws IOException {
		int nrbins = 10;
		int[][] dists = new int[files.length][nrbins];
		int[][] distsnom = new int[files.length][nrbins];
		for (int f = 0; f < files.length; f++) {
			TextFile tf = new TextFile(files[f], TextFile.R);
			String[] header = tf.readLineElems(TextFile.tab);
			
			int genecol = 0;
			int fdrcol = 0;
			int pvalcol = 0;
			int pvalthresholdcol = 0;
			for (int i = 0; i < header.length; i++) {
				String v = header[i];
				if (v.equals("gene_id")) {
					genecol = i;
				} else if (v.equals("pval_nominal")) {
					pvalcol = i;
				} else if (v.equals("qval")) {
					fdrcol = i;
				} else if (v.equals("pval_nominal_threshold")) {
					pvalthresholdcol = i;
				}
			}
			String[] elems = tf.readLineElems(TextFile.tab);
			double maxfdrpval = 0;
			
			while (elems != null) {
				
				double pval = Double.parseDouble(elems[pvalcol]);
				double log10p = -Math.log10(pval);
				double pvalnomthreshold = Double.parseDouble(elems[pvalthresholdcol]);
				
				
				if (log10p >= nrbins) {
					log10p = nrbins - 1;
				}
				
				int bin = (int) Math.floor(log10p);
				
				
				double fdr = Double.parseDouble(elems[fdrcol]);
				
				if (fdr < 0.05 && pval > maxfdrpval) {
					maxfdrpval = pval;
				}
				if (fdr < 0.05) {
					dists[f][bin]++;
				}
				if (pval < pvalnomthreshold) {
					distsnom[f][bin]++;
				}
				
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			System.out.println(files[f] + "\t" + maxfdrpval);
		}
		
		TextFile outf = new TextFile(out, TextFile.W);
		String header = "file";
		for (int i = 0; i < nrbins; i++) {
			header += "\t1E-" + i;
		}
		outf.writeln(header);
		for (int f = 0; f < files.length; f++) {
			String ln = new File(files[f]).getName() + "\t" + Strings.concat(dists[f], Strings.tab);
			outf.writeln(ln);
		}
		outf.close();
		
		TextFile outf2 = new TextFile(out + "-nominal.txt", TextFile.W);
		outf2.writeln(header);
		for (int f = 0; f < files.length; f++) {
			String ln = new File(files[f]).getName() + "\t" + Strings.concat(distsnom[f], Strings.tab);
			outf2.writeln(ln);
		}
		outf2.close();
	}
	
	public void runCheckGTExInEQTLGen(String[] files, String eqtlgenfile, String out) throws IOException, DocumentException {
		int nrbins = 10;
		
		ArrayList<HashMap<String, GTEXeQTL>> eqtls = new ArrayList<>();
		
		for (int f = 0; f < files.length; f++) {
			HashMap<String, GTEXeQTL> data = new HashMap<>();
			
			TextFile tf = new TextFile(files[f], TextFile.R);
			String[] header = tf.readLineElems(TextFile.tab);
			
			int genecol = 0;
			int fdrcol = 0;
			int pvalcol = 0;
			int refcol = 0;
			int altcol = 0;
			int slopecol = 0;
			int pvalthresholdcol = 0;
			int chrcol = 0;
			int poscol = 0;
			int topsnpcol = 0;
			// chr
			// pos
			boolean psychencode = false;
			boolean emp = false;
			String filename = new File(files[f]).getName();
			if (filename.endsWith("egenes.txt.gz")) {
				for (int i = 0; i < header.length; i++) {
					String v = header[i];
					if (v.equals("gene_id")) {
						genecol = i;
					} else if (v.equals("pval_nominal")) {
						pvalcol = i;
					} else if (v.equals("qval")) {
						fdrcol = i;
					} else if (v.equals("pval_nominal_threshold")) {
						pvalthresholdcol = i;
					} else if (v.equals("ref")) {
						refcol = i;
					} else if (v.equals("alt")) {
						altcol = i;
					} else if (v.equals("slope")) {
						slopecol = i;
					} else if (v.equals("chr")) {
						chrcol = i;
					} else if (v.equals("pos")) {
						poscol = i;
					}
				}
			} else if (filename.contains("eQTLProbes")) {
				emp = true;
				for (int i = 0; i < header.length; i++) {
					String v = header[i];
					if (v.equals("ProbeName")) {
						genecol = i;
					} else if (v.equals("PValue")) {
						pvalcol = i;
					} else if (v.equals("FDR")) {
						fdrcol = i;
					} else if (v.equals("pval_nominal_threshold")) {
						pvalthresholdcol = i;
					} else if (v.equals("SNPType")) {
						refcol = i;
					} else if (v.equals("AlleleAssessed")) {
						altcol = i;
					} else if (v.equals("OverallZScore")) {
						slopecol = i;
					} else if (v.equals("SNPChr")) {
						chrcol = i;
					} else if (v.equals("SNPChrPos")) {
						poscol = i;
					}
				}
			} else {
				// psychencode data
				psychencode = true;
				for (int i = 0; i < header.length; i++) {
					String v = header[i];
					if (v.equals("gene_id")) {
						genecol = i;
					} else if (v.equals("nominal_pval")) {
						pvalcol = i;
					} else if (v.equals("FDR")) {
						fdrcol = i;
					} else if (v.equals("pval_nominal_threshold")) {
						pvalthresholdcol = i;
					} else if (v.equals("Ref")) {
						refcol = i;
					} else if (v.equals("Alt")) {
						altcol = i;
					} else if (v.equals("regression_slope")) {
						slopecol = i;
					} else if (v.equals("SNP_chr")) {
						chrcol = i;
					} else if (v.equals("SNP_start")) {
						poscol = i;
					} else if (v.equals("top_SNP")) {
						topsnpcol = i;
					}
				}
			}
			
			
			String[] elems = tf.readLineElems(TextFile.tab);
			double maxfdrpval = 0;
			
			
			while (elems != null) {
				
				double pval = Double.parseDouble(elems[pvalcol]);
//				double pvalnomthreshold = Double.parseDouble(elems[pvalthresholdcol]);
				double fdr = Double.parseDouble(elems[fdrcol]);
				
				if (fdr < 0.05) {
					
					String gene = elems[genecol].split("\\.")[0];
					
					
					GTEXeQTL e = new GTEXeQTL();
					if (!emp) {
						e.allele = elems[refcol] + "/" + elems[altcol];
						e.alleleAssessed = elems[altcol];
						e.z = -ZScores.pToZ(pval);
						Double slp = Double.parseDouble(elems[slopecol]);
						if (slp < 0) {
							e.z *= -1;
						}
						
						if (psychencode) {
							// check whether this is the top snp
							
							int topsnp = Integer.parseInt(elems[topsnpcol]);
							if (topsnp > 0) {
								String combo = elems[chrcol].replace("chr", "") + ":" + elems[poscol] + "-" + gene;
								data.put(combo, e);
							}
							
						} else {
							String combo = elems[chrcol] + ":" + elems[poscol] + "-" + gene;
							data.put(combo, e);
						}
					} else {
						e.allele = elems[refcol];
						e.alleleAssessed = elems[altcol];
						e.z = Double.parseDouble(elems[slopecol]);
						String combo = elems[chrcol] + ":" + elems[poscol] + "-" + gene;
						data.put(combo, e);
					}
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			System.out.println(files[f] + "\t" + data.size());
			
			
			eqtls.add(data);
		}
		
		
		TextFile tf = new TextFile(eqtlgenfile, TextFile.R);
		int[] shared = new int[files.length];
		int[] significant = new int[files.length];
		int[] samedir = new int[files.length];
		
		tf.readLineElems(TextFile.tab);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		
		int nrcols = 5;
		int nrrows = (int) Math.ceil(files.length / nrcols) + 1;
		if (nrrows == 0) {
			nrrows = 1;
		}
		Grid g = new Grid(500, 500, nrrows, nrcols, 100, 250);
		
		ScatterplotPanel[] panels = new ScatterplotPanel[files.length];
		for (int i = 0; i < files.length; i++) {
			panels[i] = new ScatterplotPanel(1, 1);
			panels[i].setDataRange(new Range(-200, -50, 200, 50));
			panels[i].setAlpha(0.5f);
			panels[i].setPlotElems(true, false);
			
			
		}
		String ln = tf.readLine();
		while (ln != null) {
			elems = Strings.tab.split(ln);
			String eqtl = elems[2] + ":" + elems[3] + "-" + elems[4];
			for (int i = 0; i < files.length; i++) {
				
				HashMap<String, GTEXeQTL> l = eqtls.get(i);
				
				GTEXeQTL e = l.get(eqtl);
				Double fdr = Double.parseDouble(elems[elems.length - 1]);
				
				if (e != null) {
					shared[i]++;
					if (fdr < 0.05) {
						significant[i]++;
						String alleles = elems[8];
						String assessed = elems[9];
						Double origz = Double.parseDouble(elems[10]);
						
						Boolean flip = BaseAnnot.flipalleles(e.allele, e.alleleAssessed, alleles, assessed);
						if (flip != null) {
							double z = origz;
							if (flip) {
								z *= -1;
							}
							if (z * e.z >= 0) {
								samedir[i]++;
							}
							panels[i].addData(z, e.z);
						}
						
						
					}
				}
				
			}
			
			ctr++;
			if (ctr % 10000000 == 0) {
				System.out.println(ctr + " lines processed");
			}
			ln = tf.readLine();
		}
		tf.close();
		
		
		TextFile tfo = new TextFile(out, TextFile.W);
		tfo.writeln("File\tTotal\tShared\tSignificant\tSameDir");
		for (int i = 0; i < files.length; i++) {
			
			String name = new File(files[i]).getName();
			name = name.replace(".v7.egenes.txt.gz", "");
			name = name.replace("-DER-08a_hg19_eQTL.significant-withalleles.txt", "");
			if (name.contains("2018-01-31-cis-eQTLProbes")) {
				name = "eQTLGen";
			} else if (name.equals("eQTLProbesFDR0.05-ProbeLevel.txt.gz")) {
				name = "MetaBrain";
			}
			ln = name + "\t" + eqtls.get(i).size() + "\t" + shared[i] + "\t" + significant[i] + "\t" + samedir[i];
			
			
			double percShared = 0d;
			if (samedir[i] > 0 && significant[i] > 0) {
				percShared = ((double) samedir[i] / significant[i]) * 100;
			}
			DecimalFormat format = new DecimalFormat("##.##");
			panels[i].setTitle(eqtls.get(i).size() + " total, " + shared[i] + " shared, " + significant[i] + " significant, " + samedir[i] + " same direction (" + format.format(percShared) + "%)");
			panels[i].setLabels("eQTLGen", name);
			tfo.writeln(ln);
			
			g.addPanel(panels[i]);
		}
		tfo.close();
		g.draw(out + "ZScore-ComparisonPlot.pdf");
	}
	
	
	class GTEXeQTL {
		Double z;
		String allele;
		String alleleAssessed;
		
	}
	
	public void tssdistance(String[] files, String eqtlgenannot, String metabrainannot, String outplot) throws IOException, DocumentException {
		ArrayList<ArrayList<TSSDistance>> alldists = new ArrayList<>();
		
		int maxnrsig = 0;
		for (int f = 0; f < files.length; f++) {
			
			TextFile tf = new TextFile(files[f], TextFile.R);
			String[] header = tf.readLineElems(TextFile.tab);
			
			int genecol = 0;
			int fdrcol = 0;
			int pvalcol = 0;
			int refcol = 0;
			int altcol = 0;
			int slopecol = 0;
			int pvalthresholdcol = 0;
			int chrcol = 0;
			int poscol = 0;
			int topsnpcol = 0;
			int tssdistcol = 0;
			// chr
			// pos
			boolean psychencode = false;
			boolean emp = false;
			HashMap<String, Gene> strToGene = null;
			String filename = new File(files[f]).getName();
			if (filename.endsWith("egenes.txt.gz")) {
				for (int i = 0; i < header.length; i++) {
					String v = header[i];
					if (v.equals("gene_id")) {
						genecol = i;
					} else if (v.equals("pval_nominal")) {
						pvalcol = i;
					} else if (v.equals("qval")) {
						fdrcol = i;
					} else if (v.equals("pval_nominal_threshold")) {
						pvalthresholdcol = i;
					} else if (v.equals("ref")) {
						refcol = i;
					} else if (v.equals("alt")) {
						altcol = i;
					} else if (v.equals("slope")) {
						slopecol = i;
					} else if (v.equals("chr")) {
						chrcol = i;
					} else if (v.equals("pos")) {
						poscol = i;
					} else if (v.equals("tss_distance")) {
						tssdistcol = i;
					}
				}
			} else if (filename.contains("eQTLProbes")) {
				emp = true;
				GTFAnnotation annot;
				
				if (filename.contains("CohortInfoRemoved")) {
					annot = new GTFAnnotation(eqtlgenannot);
				} else {
					annot = new GTFAnnotation(metabrainannot);
				}
				strToGene = new HashMap<String, Gene>();
				for (Gene g : annot.getGenes()) {
					strToGene.put(g.getName(), g);
				}
				for (int i = 0; i < header.length; i++) {
					String v = header[i];
					if (v.equals("ProbeName")) {
						genecol = i;
					} else if (v.equals("PValue")) {
						pvalcol = i;
					} else if (v.equals("FDR")) {
						fdrcol = i;
					} else if (v.equals("pval_nominal_threshold")) {
						pvalthresholdcol = i;
					} else if (v.equals("SNPType")) {
						refcol = i;
					} else if (v.equals("AlleleAssessed")) {
						altcol = i;
					} else if (v.equals("OverallZScore")) {
						slopecol = i;
					} else if (v.equals("SNPChr")) {
						chrcol = i;
					} else if (v.equals("SNPChrPos")) {
						poscol = i;
					}
				}
			} else {
				// psychencode data
				psychencode = true;
				for (int i = 0; i < header.length; i++) {
					String v = header[i];
					if (v.equals("gene_id")) {
						genecol = i;
					} else if (v.equals("nominal_pval")) {
						pvalcol = i;
					} else if (v.equals("FDR")) {
						fdrcol = i;
					} else if (v.equals("pval_nominal_threshold")) {
						pvalthresholdcol = i;
					} else if (v.equals("Ref")) {
						refcol = i;
					} else if (v.equals("Alt")) {
						altcol = i;
					} else if (v.equals("regression_slope")) {
						slopecol = i;
					} else if (v.equals("SNP_chr")) {
						chrcol = i;
					} else if (v.equals("SNP_start")) {
						poscol = i;
					} else if (v.equals("top_SNP")) {
						topsnpcol = i;
					} else if (v.equals("SNP_distance_to_TSS")) {
						tssdistcol = i;
					}
				}
			}
			
			TextFile tfo = null;
			if (emp) {
				tfo = new TextFile(outplot + "-" + filename + ".txt", TextFile.W);
				
				
			}
			
			
			String[] elems = tf.readLineElems(TextFile.tab);
			double maxfdrpval = 0;
			
			ArrayList<TSSDistance> dists = new ArrayList<>();
			
			while (elems != null) {
				
				double pval = Double.parseDouble(elems[pvalcol]);
//				double pvalnomthreshold = Double.parseDouble(elems[pvalthresholdcol]);
				double fdr = Double.parseDouble(elems[fdrcol]);
				
				if (fdr < 0.05) {
					
					if (emp) {
						
						String gene = elems[genecol];
						String snp = elems[1];
						Gene geneobj = strToGene.get(gene);
						if (geneobj != null) {
							int tss = geneobj.getStart();
							if (geneobj.getStrand().equals(Strand.NEG)) {
								tss = geneobj.getStop();
							}
							int snppos = Integer.parseInt(elems[poscol]);
							int distance = snppos - tss;
							TSSDistance d = new TSSDistance();
							d.distance = distance;
							d.z = Math.abs(Double.parseDouble(elems[slopecol]));
							dists.add(d);
							tfo.writeln(pval + "\t" + gene + "\t" + snp + "\t" + distance);
						} else {
							System.out.println("Error: " + gene + " not found in annotation?");
						}
					} else {
						
						double z = -ZScores.pToZ(pval);
						if (psychencode) {
							// check whether this is the top snp
							
							int topsnp = Integer.parseInt(elems[topsnpcol]);
							if (topsnp > 0) {
								TSSDistance d = new TSSDistance();
								int dist = Integer.parseInt(elems[tssdistcol]);
								d.distance = dist;
								d.z = z;
								dists.add(d);
							}
							
						} else {
							TSSDistance d = new TSSDistance();
							int dist = Integer.parseInt(elems[tssdistcol]);
							d.distance = dist;
							d.z = z;
							dists.add(d);
						}
					}
					
					
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			if (dists.size() > maxnrsig) {
				maxnrsig = dists.size();
			}
			
			Collections.sort(dists);
			alldists.add(dists);
			if (emp) {
				tfo.close();
			}
		}
		
		
		Grid g = new Grid(500, 200, (files.length / 2) + 1, 3, 100, 100);
		for (int i = 0; i < files.length; i++) {
			
			
			ArrayList<TSSDistance> dists = alldists.get(i);
			ScatterplotPanel p = new ScatterplotPanel(1, 1);
//			p.setDataRange(new Range(0, -1000000, dists.size(), 1000000));
			p.setDataRange(new Range(0, -1000000, 35000, 1000000));
			p.setPlotElems(true, false);
			for (int d = 0; d < dists.size(); d++) {
				p.addData(d, dists.get(d).distance);
			}
			
			String name = new File(files[i]).getName();
			name = name.replace(".v7.egenes.txt.gz", "");
			name = name.replace("-DER-08a_hg19_eQTL.significant-withalleles.txt", "");
			if (name.contains("2018-01-31-cis-eQTLProbes")) {
				name = "eQTLGen";
			} else if (name.equals("eQTLProbesFDR0.05-ProbeLevel.txt.gz")) {
				name = "MetaBrain";
			}
			
			p.setTitle(name + " (" + dists.size() + " significant eQTLs)");
			p.setAlpha(0.2f);
			p.setLabels("eQTL Significance", "TSS Distance");
			
			g.addPanel(p);
		}
		g.draw(outplot);
		
		
	}
	
	public void tssdistancePerFDRBin(String[] files, String eqtlgenannot, String metabrainannot, String outplot) throws IOException, DocumentException {
		ArrayList<ArrayList<TSSDistance>> alldists = new ArrayList<>();
		int nrfdrbins = 1;
		int nrdistbins = 10;
		double[][][] fdrbins = new double[nrfdrbins][files.length][nrdistbins + 1];
		for (int f = 0; f < files.length; f++) {
			
			TextFile tf = new TextFile(files[f], TextFile.R);
			String[] header = tf.readLineElems(TextFile.tab);
			
			int genecol = 0;
			int fdrcol = 0;
			int pvalcol = 0;
			int refcol = 0;
			int altcol = 0;
			int slopecol = 0;
			int pvalthresholdcol = 0;
			int chrcol = 0;
			int poscol = 0;
			int topsnpcol = 0;
			int tssdistcol = 0;
			// chr
			// pos
			boolean psychencode = false;
			boolean emp = false;
			HashMap<String, Gene> strToGene = null;
			String filename = new File(files[f]).getName();
			if (filename.endsWith("egenes.txt.gz")) {
				for (int i = 0; i < header.length; i++) {
					String v = header[i];
					if (v.equals("gene_id")) {
						genecol = i;
					} else if (v.equals("pval_nominal")) {
						pvalcol = i;
					} else if (v.equals("qval")) {
						fdrcol = i;
					} else if (v.equals("pval_nominal_threshold")) {
						pvalthresholdcol = i;
					} else if (v.equals("ref")) {
						refcol = i;
					} else if (v.equals("alt")) {
						altcol = i;
					} else if (v.equals("slope")) {
						slopecol = i;
					} else if (v.equals("chr")) {
						chrcol = i;
					} else if (v.equals("pos")) {
						poscol = i;
					} else if (v.equals("tss_distance")) {
						tssdistcol = i;
					}
				}
			} else if (filename.contains("eQTLProbes")) {
				emp = true;
				GTFAnnotation annot;
				
				if (filename.contains("CohortInfoRemoved")) {
					annot = new GTFAnnotation(eqtlgenannot);
				} else {
					annot = new GTFAnnotation(metabrainannot);
				}
				strToGene = new HashMap<String, Gene>();
				for (Gene g : annot.getGenes()) {
					strToGene.put(g.getName(), g);
				}
				for (int i = 0; i < header.length; i++) {
					String v = header[i];
					if (v.equals("ProbeName")) {
						genecol = i;
					} else if (v.equals("PValue")) {
						pvalcol = i;
					} else if (v.equals("FDR")) {
						fdrcol = i;
					} else if (v.equals("pval_nominal_threshold")) {
						pvalthresholdcol = i;
					} else if (v.equals("SNPType")) {
						refcol = i;
					} else if (v.equals("AlleleAssessed")) {
						altcol = i;
					} else if (v.equals("OverallZScore")) {
						slopecol = i;
					} else if (v.equals("SNPChr")) {
						chrcol = i;
					} else if (v.equals("SNPChrPos")) {
						poscol = i;
					}
				}
			} else {
				// psychencode data
				psychencode = true;
				for (int i = 0; i < header.length; i++) {
					String v = header[i];
					if (v.equals("gene_id")) {
						genecol = i;
					} else if (v.equals("nominal_pval")) {
						pvalcol = i;
					} else if (v.equals("FDR")) {
						fdrcol = i;
					} else if (v.equals("pval_nominal_threshold")) {
						pvalthresholdcol = i;
					} else if (v.equals("Ref")) {
						refcol = i;
					} else if (v.equals("Alt")) {
						altcol = i;
					} else if (v.equals("regression_slope")) {
						slopecol = i;
					} else if (v.equals("SNP_chr")) {
						chrcol = i;
					} else if (v.equals("SNP_start")) {
						poscol = i;
					} else if (v.equals("top_SNP")) {
						topsnpcol = i;
					} else if (v.equals("SNP_distance_to_TSS")) {
						tssdistcol = i;
					}
				}
			}
			
			
			String[] elems = tf.readLineElems(TextFile.tab);
			double maxfdrpval = 0;
			
			while (elems != null) {
				
				double pval = Double.parseDouble(elems[pvalcol]);
//				double pvalnomthreshold = Double.parseDouble(elems[pvalthresholdcol]);
				double fdr = Double.parseDouble(elems[fdrcol]);
				
				if (fdr < 0.05) {
					
					if (emp) {
						
						String gene = elems[genecol];
						String snp = elems[1];
						Gene geneobj = strToGene.get(gene);
						if (geneobj != null) {
							int tss = geneobj.getStart();
							if (geneobj.getStrand().equals(Strand.NEG)) {
								tss = geneobj.getStop();
							}
							int snppos = Integer.parseInt(elems[poscol]);
							int distance = Math.abs(snppos - tss);
							
							int fdrbin = 0; // (int) Math.floor(((fdr / 0.05) * nrfdrbins));
							System.out.println(fdr + "\t" + fdrbin);
							if (fdrbin < 0) {
								fdrbin = 0;
							}
							if (fdrbin >= nrfdrbins) {
								fdrbin = nrfdrbins - 1;
							}
							if (distance > 500000) {
							
							}
//							if (distance < -500000) {
//								distance = -500000;
//							}
							// distance = distance; // set -500k to 0
							double distperc = (double) distance / 500000;
							int distbin = (int) Math.floor(distperc * nrdistbins);
							if (distance > 500000) {
								distbin = nrdistbins;
							}
							
							if (distbin < 0) {
								distbin = 0;
							}
							fdrbins[fdrbin][f][distbin]++;
						} else {
							System.out.println("Error: " + gene + " not found in annotation?");
						}
					} else {
						boolean istopsnp = false;
						
						double z = -ZScores.pToZ(pval);
						int dist = 0;
						if (psychencode) {
							// check whether this is the top snp
							
							int topsnp = Integer.parseInt(elems[topsnpcol]);
							if (topsnp > 0) {
								istopsnp = true;
								dist = Integer.parseInt(elems[tssdistcol]);
								
							}
							
						} else {
							
							dist = Integer.parseInt(elems[tssdistcol]);
							istopsnp = true;
						}
						
						if (istopsnp) {
							int snppos = Integer.parseInt(elems[poscol]);
							int distance = Math.abs(dist);
							
							int fdrbin = 0; // (int) Math.floor(((fdr / 0.05) * nrfdrbins));
							if (fdrbin < 0) {
								fdrbin = 0;
							}
							if (fdrbin >= nrfdrbins) {
								fdrbin = nrfdrbins - 1;
							}
							if (distance > 500000) {
								distance = 500000;
							}
//							if (distance < -500000) {
//								distance = -500000;
//							}
							// distance = distance + 500000; // set -500k to 0
							double distperc = (double) distance / 500000;
							int distbin = (int) Math.floor(distperc * nrdistbins);
							if (distance > 500000) {
								distbin = nrdistbins;
							}
							
							if (distbin < 0) {
								distbin = 0;
							}
							fdrbins[fdrbin][f][distbin]++;
						}
					}
					
					
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			
		}
		
		
		Grid g = new Grid(500, 500, nrfdrbins, 1, 100, 100);
		
		for (int i = 0; i < nrfdrbins; i++) {
			for (int d = 0; d < fdrbins[i].length; d++) {
				double sum = 0;
				for (int j = 0; j < fdrbins[i][d].length; j++) {
					sum += fdrbins[i][d][j];
				}
				for (int j = 0; j < fdrbins[i][d].length; j++) {
					fdrbins[i][d][j] = fdrbins[i][d][j] / sum;
				}
				
			}
		}
		
		
		String[] binlabels = new String[nrdistbins+1];
		for (int i = 0; i < nrdistbins+1; i++) {
			int pos = (500000 / nrdistbins) * i;
			binlabels[i] = "" + pos;
		}
		
		for (int i = 0; i < nrfdrbins; i++) {
			HistogramPanel p = new HistogramPanel(1, 1);
			p.setData(fdrbins[i]);
			
			p.setAxisLabels("Distance to TSS", "Proportion of eQTLs");
			p.setRange(new Range(0d, 0d, 500000d, 1d));
			double thisinc = ((0.05 / nrfdrbins) * i);
			double previnc = ((0.05 / nrfdrbins) * (i - 1));
			if (i == 0) {
				p.setTitle("FDR = 0");
			} else {
				p.setTitle(previnc + " < FDR < " + thisinc);
			}
			p.setBinLabels(binlabels);
			g.addPanel(p);
		}
		
		g.draw(outplot);
		
		
	}
	
	class TSSDistance implements Comparable<TSSDistance> {
		int distance;
		double z;
		
		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			
			TSSDistance that = (TSSDistance) o;
			
			return Double.compare(that.z, z) == 0;
		}
		
		@Override
		public int hashCode() {
			long temp = Double.doubleToLongBits(z);
			return (int) (temp ^ (temp >>> 32));
		}
		
		@Override
		public int compareTo(TSSDistance o) {
			if (o.equals(this)) {
				return 0;
			} else if (o.z > z) {
				return 1;
			} else {
				return -1;
			}
		}
	}
	
	
}
