package nl.harmjanwestra.playground.conv;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;

import java.io.*;
import java.util.ArrayList;

public class ConvertVCFToTT {
	
	
	public static void main(String[] args) {
		ConvertVCFToTT v = new ConvertVCFToTT();
		boolean writegenotypes = true;
		boolean parsevcf = true;
		try {
			v.run(args[0], args[1], parsevcf, writegenotypes);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void run(String vcffile, String output, boolean parseVCF, boolean writegenotypes) throws IOException {
		
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
		
		
		TextFile tf = new TextFile(vcffile, TextFile.R);
		String ln = tf.readLine();
		int ctr = 0;
		int written = 0;
		int nrvarswithdosage = 0;
		while (ln != null) {
			if (!ln.startsWith("#")) {
				
				VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.ALL);
				
				if (!var.isIndel() && var.getAlleles().length == 2 && var.getMAF() > 0.001) {
					
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
					
					int missing = 0;
					for (int i = 0; i < samples.size(); i++) {
						double a1 = genotypes.getQuick(i, 0);
						double a2 = genotypes.getQuick(i, 1);
						
						if (a1 == 0 || a2 == 0 || (apdepth != null && apdepth[i] < 10)) {
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
					
					if ((double) missing / samples.size() > 0.05) {
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
					}
					
					
				}
				
			}
			if (ctr % 1000 == 0) {
				System.out.print(ctr + "\tlines read. " + written + "\twritten. " + nrvarswithdosage + "\twith dosages\r");
				
			}
			ctr++;
			ln = tf.readLine();
		}
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
