package nl.harmjanwestra.playground.conv;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.*;
import java.util.ArrayList;

public class ConvertVCFToTT {
	
	
	public static void main(String[] args) {
		ConvertVCFToTT v = new ConvertVCFToTT();
		boolean writegenotypes = true;
		boolean parsevcf = true;
		boolean useAD = false;
		try {
			v.run(args[0], args[1], parsevcf, writegenotypes, useAD);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void run(String vcffile, String output, boolean parseVCF, boolean writegenotypes, boolean useAD) throws IOException {
		
		VCFGenotypeData d = new VCFGenotypeData(vcffile);
		ArrayList<String> samples = d.getSamples();
		d.close();
		
		
		String folder = output;
		File genotypeDataFile = new File(folder, "GenotypeMatrix.dat");
		File imputedDosageDataFile = new File(folder, "ImputedDosageMatrix.dat");
		File snpFile = new File(folder, "SNPs.txt.gz");
		File snpMapFile = new File(folder, "SNPMappings.txt.gz");
		File individualFile = new File(folder, "Individuals.txt");
		File phenotypeAnnotationFile = new File(folder, "PhenotypeInformation.txt");
		File logfile = new File(folder, "ConversionLog.txt.gz");
		
		BufferedWriter phenowriter = new BufferedWriter(new FileWriter(phenotypeAnnotationFile));
		TextFile indsout = new TextFile(individualFile, TextFile.W);
		for (int s = 0; s < samples.size(); s++) {
			indsout.writeln(samples.get(s));
			phenowriter.write(samples.get(s) + "\tunknown\tinclude\tunknown\n");
		}
		indsout.close();
		phenowriter.close();
		if (!parseVCF) {
			System.exit(0);
		}
		
		BufferedOutputStream genotypeDataFileWriter = null;
		BufferedOutputStream genotypeDosageDataFileWriter = null;
		
		if (writegenotypes) {
			genotypeDataFileWriter = new BufferedOutputStream(new FileOutputStream(genotypeDataFile), 32 * 1024);
			genotypeDosageDataFileWriter = new BufferedOutputStream(new FileOutputStream(imputedDosageDataFile), 32 * 1024);
		}
		TextFile snpFileWriter = new TextFile(snpFile, TextFile.W);
		TextFile snpMapFileWriter = new TextFile(snpMapFile, TextFile.W);
		
		TextFile logout = new TextFile(logfile, TextFile.W);
		
		TextFile tf = new TextFile(vcffile, TextFile.R);
		String ln = tf.readLine();
		int ctr = 0;
		int written = 0;
		int variantswithmissingness = 0;
		int nrvarswithdosage = 0;
		int nrindels = 0;
		int nrmultiallelic = 0;
		int lowmaf = 0;
		double mafthresold = 0.001;
		double missingnessthreshold = 0.05;
		int allelicdepththreshold = 10;
		String header = "var\talleles\tisIndel\tisMultiAllelic\tMAF\tMissingness\tnrMissing\tAvgDepth\tVarIncluded";
		logout.writeln(header);
		while (ln != null) {
			if (!ln.startsWith("#")) {
				
				VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.ALL);
				double missingness = 0;
				int missing = 0;
				
				int nrcalled = 0;
				boolean varincluded = false;
				double sumdepth = 0;
				if (!var.isIndel() && var.getAlleles().length == 2 && var.getMAF() > mafthresold) {
					
					// if it has read depth, set <10 as missing
					short[] apdepth = var.getApproximateDepth();
					byte allele1 = BaseAnnot.toByte(var.getAlleles()[0]);
					byte allele2 = BaseAnnot.toByte(var.getAlleles()[1]);
					
					
					DoubleMatrix2D genotypes = var.getGenotypeAllelesAsMatrix2D();
					DoubleMatrix2D dosages = var.getDosagesAsMatrix2D();
					String id = var.getId();
					if (id.length() == 1) {
						id = var.getChr() + ":" + var.getPos();
					}
					
					
					byte[] al1 = new byte[samples.size() * 2];
					
					byte[] dos = new byte[samples.size()];
					if (var.hasImputationDosages()) {
						nrvarswithdosage++;
					}
					
					missing = 0;
					for (int i = 0; i < samples.size(); i++) {
						double a1 = genotypes.getQuick(i, 0);
						double a2 = genotypes.getQuick(i, 1);
						if (apdepth != null) {
							sumdepth += apdepth[i];
						}
						if (a1 == -1 || a2 == -1 || (useAD && (apdepth != null && apdepth[i] < allelicdepththreshold))) {
							missing++;
						} else {
							if (a1 == 0) {
								al1[i] = allele1;
							} else if (a1 == 1) {
								al1[i] = allele2;
							}
							if (a2 == 0) {
								al1[i + samples.size()] = allele1;
							} else if (a2 == 1) {
								al1[i + samples.size()] = allele2;
							}
							if (var.hasImputationDosages()) {
								dos[i] = (byte) dosages.getQuick(i, 0);
							}
						}
					}
					sumdepth /= samples.size();
					missingness = (double) missing / samples.size();
					if (missingness < missingnessthreshold) {
						snpFileWriter.append(id);
						snpFileWriter.append('\n');
						snpMapFileWriter.append(var.getChr());
						snpMapFileWriter.append('\t');
						snpMapFileWriter.append(String.valueOf(var.getPos()));
						snpMapFileWriter.append('\t');
						snpMapFileWriter.append(id);
						snpMapFileWriter.append('\n');
						
						if (writegenotypes) {
							genotypeDataFileWriter.write(al1);
							genotypeDosageDataFileWriter.write(dos);
						}
						written++;
						varincluded = true;
					} else {
						variantswithmissingness++;
					}
				} else if (var.getMAF() < mafthresold) {
					lowmaf++;
				} else if (var.isIndel()) {
					nrindels++;
				} else if (var.isMultiallelic()) {
					nrmultiallelic++;
				}
				
				String incvar = "NotWritten";
				if (varincluded) {
					incvar = "Written";
				}
				String log = var.toString() + "\t" + Strings.concat(var.getAlleles(), Strings.semicolon) + "\t" + var.isIndel() + "\t" + var.isMultiallelic() + "\t" + var.getMAF() + "\t" + missingness + "\t" + missing + "\t" + sumdepth + "\t" + incvar;
				logout.writeln(log);
			}
			if (ctr % 1000 == 0) {
				System.out.print(ctr + "\tlines read. " + written + " written, " + nrvarswithdosage + " with dosages. " + lowmaf + " maf<" + mafthresold + ", " + nrindels + " indels ," + nrmultiallelic + " multiallelic, " + variantswithmissingness + " with MAF>" + mafthresold + " and high missigness\r");
				
			}
			ctr++;
			ln = tf.readLine();
		}
		logout.close();
		tf.close();
		System.out.println();
		System.out.println("Done");
		if (writegenotypes) {
			genotypeDataFileWriter.close();
			genotypeDosageDataFileWriter.close();
			if (nrvarswithdosage == 0) {
				Gpio.delete(folder + "ImputedDosageMatrix.dat");
			}
		}
		snpFileWriter.close();
		snpMapFileWriter.close();
		
		
	}
	
	
}
