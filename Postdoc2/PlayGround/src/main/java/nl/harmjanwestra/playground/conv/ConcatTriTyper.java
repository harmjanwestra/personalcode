package nl.harmjanwestra.playground.conv;

import nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

public class ConcatTriTyper {
	
	public static void main(String[] args) {
		if (args.length < 2) {
			System.out.println("Usage: concat.jar ./inputCHR/ ./output/");
		} else {
			ConcatTriTyper tt = new ConcatTriTyper();
			try {
				tt.run(args[0], args[1]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public void run(String input, String output) throws IOException {
		Gpio.createDir(output);
		BinaryFile bf1 = new BinaryFile(output + "GenotypeMatrix.dat", BinaryFile.W);
		BinaryFile bfimp1 = new BinaryFile(output + "ImputedDosageMatrix.dat", BinaryFile.W);
		BinaryFile bf2 = new BinaryFile(output + "GenotypeMatrixV2.dat", BinaryFile.W);
		BinaryFile bfimp2 = new BinaryFile(output + "ImputedDosageMatrixV2.dat", BinaryFile.W);
		TextFile snpo1 = new TextFile(output + "SNPs.txt.gz", TextFile.W);
		TextFile snpo2 = new TextFile(output + "SNPsV2.txt.gz", TextFile.W);
		TextFile snpmap = new TextFile(output + "SNPMappings.txt.gz", TextFile.W);
		TextFile convlog = new TextFile(output + "ConversionLog.txt.gz", TextFile.W);
		
		for (int i = 1; i < 26; i++) {
			
			
			String inputChr = null;
			if (i < 23) {
				inputChr = input.replaceAll("CHR", "" + i);
			} else {
				if (i == 23) {
					inputChr = input.replaceAll("CHR", "X");
				} else if (i == 24) {
					inputChr = input.replaceAll("CHR", "Y");
				} else if (i == 25) {
					inputChr = input.replaceAll("CHR", "N");
				}
				
			}
			
			String input1 = inputChr + "GenotypeMatrix.dat";
			if (!Gpio.exists(input1)) {
				System.out.println("Skipping chr: " + i);
			} else {
				String indfile = inputChr + "Individuals.txt";
				Gpio.copyFile(indfile, output + "Individuals.txt");
				Gpio.copyFile(inputChr + "PhenotypeInformation.txt", output + "PhenotypeInformation.txt");
				
				concatText(inputChr + "SNPMappings.txt.gz", snpmap);
				concatText(inputChr + "ConversionLog.txt.gz", convlog);
				
				
				TextFile ind = new TextFile(indfile, TextFile.R);
				String[] individuals = ind.readAsArray();
				
				int nrsnps = concatText(inputChr + "SNPs.txt.gz", snpo1);
				
				concatBinary(input1, bf1, individuals.length, nrsnps, false);
				String imp1 = inputChr + "ImputedDosageMatrix.dat";
				if (Gpio.exists(imp1)) {
					concatBinary(imp1, bfimp1, individuals.length, nrsnps, true);
				}
				
				
				String input2 = inputChr + "GenotypeMatrixV2.dat";
				if (Gpio.exists(input2)) {
					int nrsnps2 = concatText(inputChr + "SNPsV2.txt.gz", snpo2);
					
					concatBinary(input2, bf2, individuals.length, nrsnps2, false);
					String imp2 = inputChr + "ImputedDosageMatrixV2.dat";
					if (Gpio.exists(imp2)) {
						concatBinary(imp2, bfimp2, individuals.length, nrsnps2, true);
					}
				}
			}
			
			
		}
		
		
		bf1.close();
		bf2.close();
		bfimp1.close();
		bfimp2.close();
		snpo1.close();
		snpo2.close();
		snpmap.close();
		convlog.close();
		
		
	}
	
	private int concatText(String s, TextFile snpmap) throws IOException {
		System.out.println("Merging " + s + " into " + snpmap.getFileName());
		if (Gpio.exists(s)) {
			TextFile tc = new TextFile(s, TextFile.R);
			String snpmapln = tc.readLine();
			int nrlns = 0;
			while (snpmapln != null) {
				snpmap.writeln(snpmapln);
				snpmapln = tc.readLine();
				nrlns++;
			}
			tc.close();
			return nrlns;
		} else {
			System.out.println("Warning: could not find file: " + s);
			System.exit(-1);
			return 0;
		}
	}
	
	private void concatBinary(String input1, BinaryFile bf1, int nrbytesperAllele, int nrsnps, boolean imputed) throws IOException {
		BinaryFile bf = new BinaryFile(input1, BinaryFile.R);
		byte[] buffer = new byte[nrbytesperAllele * 2];
		if (imputed) {
			buffer = new byte[nrbytesperAllele];
		}
		
		ProgressBar pb = new ProgressBar(nrsnps, "Importing " + input1);
		for (int i = 0; i < nrsnps; i++) {
			int nrbytes = bf.read(buffer);
			bf1.write(buffer, 0, nrbytes);
			pb.iterate();
		}
		pb.close();
		
		
	}
}
