package nl.harmjanwestra.playground.cis;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.*;

public class LDFromZMat {
	
	
	public static void main(String[] args) {
		LDFromZMat z = new LDFromZMat();
		try {
//			z.createCombos("D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\sortedGeneSNPCombos.txt.gz",
//					"D:\\snpgenecombos\\");
			
			if (args.length >= 7) {
				String genecovariancematrixloc = args[0];
				String cohortselectionfile = args[1];
				String zmatloc = args[2];
				String genedefloc = args[3];
				String genelistfile = args[4];
				String referenceallelefile = args[5];
				String outloc = args[6];
				Integer minnrvals = null;
				if (args.length >= 8) {
					minnrvals = Integer.parseInt(args[7]);
				}
				z.calculateLD(genecovariancematrixloc, cohortselectionfile, zmatloc, genedefloc, genelistfile, referenceallelefile, outloc, minnrvals);
			} else {
				System.out.println("Arguments:\ngenecovariancematrix cohortselectionfile zmatloc genedefloc genelistfile referenceallelefile outploc [minnrvals]");
				
			}
			
			System.exit();
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void createCombos(String genesnpcombos, String snpgeneout) throws IOException {
		
		
		TextFile tf = new TextFile(genesnpcombos, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		
		
		LinkedHashMap<String, Short> geneMap = new LinkedHashMap<String, Short>();
		
		HashMap<String, ArrayList<Short>> snpgenecombos = new HashMap<>();
		
		
		short gctr = 0;
		int lctr = 0;
		while (elems != null) {
			
			String gene = Strings.cache(elems[0]);
			String snp = Strings.cache(elems[1]);
			
			if (!geneMap.containsKey(gene)) {
				geneMap.put(gene, gctr);
				gctr++;
				if (gctr < 0) {
					System.exit(-1);
				}
			}
			
			short geneId = geneMap.get(gene);
			
			
			ArrayList<Short> genes = snpgenecombos.get(snp);
			if (genes == null) {
				genes = new ArrayList<>();
			}
			genes.add(geneId);
			snpgenecombos.put(snp, genes);
			
			lctr++;
			if (lctr % 100000 == 0) {
				System.out.print("\r" + lctr + " lines parsed ");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		
		System.out.println(snpgenecombos.size() + " snps loaded. ");
		tf.open();
		
		
		ArrayList<String> geneList = new ArrayList<>();
		geneList.addAll(geneMap.keySet());
		
		elems = tf.readLineElems(TextFile.tab);
		String prevgene = null;
		TextFile out = null;
		lctr = 0;
		while (elems != null) {
			String gene = Strings.cache(elems[0]);
			String snp = Strings.cache(elems[1]);
			
			if (prevgene == null || !gene.equals(prevgene)) {
				if (out != null) {
					out.close();
				}
				out = new TextFile(snpgeneout + gene + ".txt.gz", TextFile.W);
			}
			
			prevgene = gene;
			
			ArrayList<Short> combos = snpgenecombos.get(snp);
			String[] genes = new String[combos.size()];
			for (int i = 0; i < genes.length; i++) {
				genes[i] = geneList.get(i);
			}
			out.writeln(snp + "\t" + Strings.concat(genes, Strings.semicolon));
			
			lctr++;
			if (lctr % 100000 == 0) {
				System.out.println("Lines parsed: " + lctr + " current gene: " + prevgene);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		
		if (out != null) {
			out.close();
		}
		
	}
	
	
	public void calculateLD(String genecovariancematrixloc,
							String cohortSelectionFile,
							String zmatloc,
							String geneDefLoc,
							String geneListFile,
							String referenceAlleleFile,
							String outloc,
							Integer minNrVals) throws Exception {
		
		boolean run = true;
		if (!Gpio.exists(genecovariancematrixloc)) {
			System.out.println("Could not find gene covariance file: " + genecovariancematrixloc);
			run = false;
		}
		if (!Gpio.exists(cohortSelectionFile)) {
			System.out.println("Could not find cohort selection file: " + cohortSelectionFile);
			run = false;
		}
		if (!Gpio.exists(zmatloc)) {
			System.out.println("Could not find z-score matrix location: " + zmatloc);
			run = false;
		}
		if (!Gpio.exists(geneDefLoc)) {
			System.out.println("Could not find gene definition location: " + geneDefLoc);
			run = false;
		}
		if (!Gpio.exists(geneListFile)) {
			System.out.println("Could not find gene list file: " + geneListFile);
			run = false;
		}
		if (!Gpio.exists(referenceAlleleFile)) {
			System.out.println("Could not find reference allele file: " + referenceAlleleFile);
			run = false;
		}
		if (!run) {
			System.out.println("Some run conditions not met.");
			System.exit(-1);
		}
		if (!Gpio.exists(outloc)) {
			Gpio.createDir(outloc);
		}
		
		//		String[] cohortsEuropean = {"inCHIANTI", "HVH_HT12v3", "EGCUT_RNAseq", "LIFE_Adult_plus", "EGCUT_HT12v4", "Fehrmann_HT12v3", "NTR_NESDA", "LIFE_Heart_plus", "BSGS", "CODAM", "LLS_660Q", "DILGOM", "Fehrmann_H8v2", "GoNL_WGS", "NTR_GoNL", "SHIP_TREND", "CARTaGENE_freeze2", "Rotterdam_HT12v4", "KORA_F4", "Cardiology", "NTR_AFFY", "GTEx", "EGCUT_HT12v3", "DGN", "CHDWB", "CARTaGENE_freeze1", "PAN", "LLS_OmniExpr", "Sorbs", "FHS", "HVH_HT12v4", "Rotterdam_RNASeq", "LL"};
//		String[] genes = {"ENSG00000002726", "ENSG00000002933", "ENSG00000055118", "ENSG00000106565", "ENSG00000164867", "ENSG00000177590", "ENSG00000241134"};
//		String genecovariancematrixloc = "/Users/lude/Documents/Genetica/eQTLGen/ZScoreMarices-20180125//GeneCorrelation/ProbeCovariance-Perms1To10.binary";
//
//		String cohortSelectionFile = null;
//		String zmatloc = null;
//		String geneDefLoc = null;
//		String geneListFile = null;
//		String outloc;
		
		if (minNrVals == null) {
			minNrVals = 0;
		}
		
		System.out.println("Loading genes from: " + geneListFile);
		TextFile tfc = new TextFile(geneListFile, TextFile.R);
		ArrayList<String> genesToRun = tfc.readAsArrayList();
		tfc.close();
		
		System.out.println("Selected " + genesToRun.size() + " genes to calculate LD for.");
		if (genesToRun.isEmpty()) {
			System.out.println("No genes found to run");
			System.exit(-1);
		}
		
		HashSet<String> hashCohortsEuropean = new HashSet<String>();
		{
			System.out.println("Loading cohorts from file: " + cohortSelectionFile);
			TextFile tfcs = new TextFile(cohortSelectionFile, TextFile.R);
			ArrayList<String> cohortsToInclude = tfcs.readAsArrayList();
			tfcs.close();
			
			for (int c = 0; c < cohortsToInclude.size(); c++) {
				for (int perm = 1; perm <= 10; perm++) {
					hashCohortsEuropean.add(cohortsToInclude.get(c) + "-perm-" + perm);
				}
			}
		}
		
		HashMap<String, String> referenceAlleleMap = new HashMap<String, String>();
		{
			TextFile tf = new TextFile(referenceAlleleFile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				String snp = elems[0];
				String alleles = elems[1] + ";" + elems[2];
				referenceAlleleMap.put(new String(snp.getBytes("UTF-8")), Strings.cache(alleles));
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		}
		
		System.out.println("Requested: " + hashCohortsEuropean.size() + " cohorts.");
		if (hashCohortsEuropean.isEmpty()) {
			System.out.println("No cohorts specified");
			System.exit(-1);
		}
		
		// load reference alleles
		for (String querygene : genesToRun) {
			// determine which snps and genes to read..
			
			String geneDef = geneDefLoc + querygene + ".txt.gz";
			if (!Gpio.exists(geneDef)) {
				System.out.println("Could not find gene definition: " + geneDef);
				break;
			} else {
				LinkedHashSet<String> snpSet = new LinkedHashSet<String>();
				ArrayList<String> snpList = new ArrayList<>();
				LinkedHashSet<String> geneSet = new LinkedHashSet<String>();
				
				TextFile tf = new TextFile(geneDef, TextFile.R);
				String[] elems = tf.readLineElems(TextFile.tab);
				while (elems != null) {
					String snp = elems[0];
					
					String[] genes = elems[1].split(";");
					snpSet.add(new String(snp.getBytes("UTF-8")));
					
					for (String gene : genes) {
						geneSet.add(new String(gene.getBytes("UTF-8")));
					}
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
				
				snpList.addAll(snpSet);
				
				double[] weights = null;
				String[] genesWithWeights = null;
				{
					DoubleMatrixDataset<String, String> datasetCorr = DoubleMatrixDataset.loadSubsetOfBinaryDoubleData(
							genecovariancematrixloc, geneSet, geneSet);
					
					//Calculate the factorloadings:
					Jama.EigenvalueDecomposition eig = eigenValueDecomposition(datasetCorr.getMatrix().toArray());
					double[] eigenValues = eig.getRealEigenvalues();
					
					DoubleMatrix2D factorLoadings = new DenseDoubleMatrix2D(datasetCorr.rows(), datasetCorr.rows());
					for (int comp = 0; comp < datasetCorr.rows(); comp++) {
						double eigenvalue = eigenValues[eigenValues.length - 1 - comp];
						if (eigenvalue < 0) {
							eigenvalue = 0;
						}
						double sqrtEigenvalue = Math.sqrt(eigenvalue);
						double[] eigenvector = getEigenVector(eig, comp);
						
						for (int a = 0; a < datasetCorr.rows(); a++) {
							double v = sqrtEigenvalue * eigenvector[a];
							factorLoadings.setQuick(comp, a, v);
						}
					}
					
					//Calculate the weights of the individual genes, to be used for the weighed correlation, to account for co-expression between genes:
					weights = new double[datasetCorr.rows()];
					genesWithWeights = new String[datasetCorr.rows()];
					for (int p = 0; p < datasetCorr.rows(); p++) {
						double weight = 0;
						genesWithWeights[p] = datasetCorr.getRowObjects().get(p);
						for (int comp = 0; comp < datasetCorr.rows(); comp++) {
							double eigenvalue = eigenValues[eigenValues.length - 1 - comp];
							if (eigenvalue < 1) {
								eigenvalue = 1;
							}
							double factorLoadingCompP = factorLoadings.getQuick(comp, p);
							weight += factorLoadingCompP * factorLoadingCompP / eigenvalue;
						}
						
						weights[p] = weight;
						System.out.println(p + "\t" + datasetCorr.getRowObjects().get(p) + "\t" + weights[p]);
					}
				}
				
				ArrayList<DoubleMatrixDataset<String, String>> zmats = new ArrayList<>();
				
				
				// load relevant z-score matrices..
				for (int g = 0; g < genesWithWeights.length; g++) {
					String geneToLoad = zmatloc + genesWithWeights[g] + "-zmat.txt.gz";
					String allelefile = zmatloc + genesWithWeights[g] + "-referenceAlleles.txt.gz";
					
					if (Gpio.exists(geneToLoad)) {
						DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadSubsetOfTextDoubleData(zmatloc, '\t', hashCohortsEuropean, snpSet);
						zmats.add(ds);
						// flip alleles
						TextFile tfa = new TextFile(allelefile, TextFile.R);
						elems = tfa.readLineElems(TextFile.tab);
						while (elems != null) {
							
							String snp = elems[0];
							Integer snpid = ds.getHashCols().get(snp);
							if (snpid != null) {
								String alleles = elems[1];
								String assessed = elems[2];
								String referenceAlleles = referenceAlleleMap.get(snp);
								if (referenceAlleles != null) {
									String[] refAlleleElems = referenceAlleles.split(";");
									Boolean flip = BaseAnnot.flipalleles(refAlleleElems[0], refAlleleElems[1], alleles, assessed);
									if (flip == null) {
										// kill snp
										for (int r = 0; r < ds.rows(); r++) {
											ds.getMatrix().setQuick(r, snpid, Double.NaN);
										}
									} else if (flip) {
										// kill snp
										for (int r = 0; r < ds.rows(); r++) {
											ds.getMatrix().setQuick(r, snpid, ds.getMatrix().getQuick(r, snpid) * -1);
										}
									}
								} else {
									// kill snp
									for (int r = 0; r < ds.rows(); r++) {
										ds.getMatrix().setQuick(r, snpid, Double.NaN);
									}
								}
							}
							
							elems = tfa.readLineElems(TextFile.tab);
						}
						
						
						tfa.close();
					} else {
						zmats.add(null);
					}
				}
				
				int nrSNPs = snpSet.size();
				
				DoubleMatrixDataset<String, String> ldmatrix = new DoubleMatrixDataset<String, String>(nrSNPs, nrSNPs);
				ldmatrix.setRowObjects(snpList);
				ldmatrix.setColObjects(snpList);
				TextFile tfout = new TextFile(outloc + querygene + "-list.txt.gz", TextFile.W);
				
				DoubleMatrix2D matrixObj = ldmatrix.getMatrix();
				for (int snp1 = 0; snp1 < nrSNPs; snp1++) {
					String snp1Name = snpList.get(snp1);
					matrixObj.setQuick(snp1, snp1, 1d);
					
					for (int snp2 = snp1 + 1; snp2 < nrSNPs; snp2++) {
						ArrayList<Double> snp1d = new ArrayList<>();
						ArrayList<Double> snp2d = new ArrayList<>();
						ArrayList<Double> weightsd = new ArrayList<>();
						String snp2Name = snpList.get(snp2);
						
						// accumulate data over genes
						for (int g = 0; g < zmats.size(); g++) {
							DoubleMatrixDataset<String, String> geneZmat = zmats.get(g);
							if (geneZmat != null) {
								Integer snp1id = geneZmat.getHashCols().get(snp1Name);
								Integer snp2id = geneZmat.getHashCols().get(snp2Name);
								
								if (snp1id != null && snp2id != null) {
									for (int row = 0; row < geneZmat.rows(); row++) {
										double v1 = geneZmat.getElementQuick(row, snp1id);
										double v2 = geneZmat.getElementQuick(row, snp2id);
										
										if (!Double.isNaN(v1) && !Double.isNaN(v2)) {
											snp1d.add(v1);
											snp2d.add(v2);
											weightsd.add(weights[g]);
										}
									}
								}
							}
						}
						
						// calculate weighted correlation
						if (snp1d.size() >= minNrVals) {
							double[] vals1 = Primitives.toPrimitiveArr(snp1d);
							double[] vals2 = Primitives.toPrimitiveArr(snp2d);
							double[] valWeights = Primitives.toPrimitiveArr(weightsd);
							
							double corr = JSci.maths.ArrayMath.correlation(vals1, vals2);
							double weightedCorr = weightedCorrelation(vals1, vals2, valWeights);
							System.out.println(snp1 + "\t" + snp2 + "\t" + snp1Name + "\t" + snp2Name + "\t" + snp1d.size() + "\t" + corr + "\t" + weightedCorr);
							tfout.writeln(snp1 + "\t" + snp2 + "\t" + snp1Name + "\t" + snp2Name + "\t" + snp1d.size() + "\t" + corr + "\t" + weightedCorr);
							matrixObj.setQuick(snp1, snp2, weightedCorr);
							matrixObj.setQuick(snp2, snp1, weightedCorr);
						}
						
					}
				}
				tfout.close();
				ldmatrix.save(outloc + querygene + "-ldmatrix.txt.gz");
			}
		}
		
		
	}
	
	public double weightedCorrelation(double[] x, double[] y, double[] weights) {
		double wmX = weightedMean(x, weights);
		double wmY = weightedMean(y, weights);
		double sumWeights = JSci.maths.ArrayMath.mass(weights);
		double covXX = 0;
		double covXY = 0;
		double covYY = 0;
		for (int s = 0; s < x.length; s++) {
			covXX += weights[s] * (x[s] - wmX) * (x[s] - wmX);
			covXY += weights[s] * (x[s] - wmX) * (y[s] - wmY);
			covYY += weights[s] * (y[s] - wmY) * (y[s] - wmY);
		}
		covXX /= sumWeights;
		covXY /= sumWeights;
		covYY /= sumWeights;
		double corr = covXY / (Math.sqrt(covXX * covYY));
		return corr;
	}
	
	public double weightedMean(double[] x, double[] weights) {
		double m = 0;
		double sumWeights = 0;
		for (int s = 0; s < x.length; s++) {
			m += x[s] * weights[s];
			sumWeights += weights[s];
		}
		return m / sumWeights;
	}
	
	
	private Jama.EigenvalueDecomposition eigenValueDecomposition(double[][] data) {
		Jama.Matrix m = new Jama.Matrix(data);
		Jama.EigenvalueDecomposition eig = m.eig();
		return eig;
	}
	
	private double[] getEigenVector(Jama.EigenvalueDecomposition eig, double[] eigenValues, int pca) {
		Jama.Matrix eigenValueMatrix = eig.getV();
		double[][] eigenValueMat = eigenValueMatrix.getArray();
		double[] eigenVector = new double[eigenValueMat.length];
		for (int i = 0; i < eigenValueMat.length; i++) {
			eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
		}
		return eigenVector;
	}
	
	private double[] getEigenVector(Jama.EigenvalueDecomposition eig, int pca) {
		Jama.Matrix eigenValueMatrix = eig.getV();
		double[][] eigenValueMat = eigenValueMatrix.getArray();
		double[] eigenVector = new double[eigenValueMat.length];
		for (int i = 0; i < eigenValueMat.length; i++) {
			eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
		}
		return eigenVector;
	}
	
	
	private double getEigenValueVar(double[] eigenValues, int pca) {
		double sumEigenvalues = 0.0;
		for (Double d : eigenValues) {
			sumEigenvalues += Math.abs(d);
		}
		double result = eigenValues[eigenValues.length - 1 - pca] / sumEigenvalues;
		return result;
	}
	
	private double[] getEigenVectorSVD(Jama.SingularValueDecomposition svd, double[] singularValues, int pca) {
		Jama.Matrix eigenValueMatrix = svd.getV();
		double[][] eigenValueMat = eigenValueMatrix.getArray();
		double[] eigenVector = new double[eigenValueMat.length];
		for (int i = 0; i < eigenValueMat.length; i++) {
			eigenVector[i] = eigenValueMat[i][pca] * Math.sqrt(singularValues[pca]);
		}
		return eigenVector;
	}
}
