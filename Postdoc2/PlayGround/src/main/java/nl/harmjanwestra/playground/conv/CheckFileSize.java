package nl.harmjanwestra.playground.conv;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

public class CheckFileSize {
	
	public static void main(String[] args) {
		CheckFileSize f = new CheckFileSize();
		try {
			f.run(args[0]);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String indir) throws IOException {
		TextFile tf = new TextFile(indir + "SNPs.txt.gz", TextFile.R);
		int nrsnps = tf.countLines();
		tf.close();
		
		tf = new TextFile(indir + "Individuals.txt", TextFile.R);
		int nrinds = tf.countLines();
		tf.close();
		
		
		if (Gpio.exists(indir + "GenotypeMatrix.dat")) {
			long expsize = ((long) nrinds * 2) * (long) nrsnps;
			long actsize = Gpio.getFileSize(indir + "GenotypeMatrix.dat");
			if (expsize != actsize) {
				System.out.println(indir + "GenotypeMatrix.dat: " + nrsnps + " snps, " + nrinds + " individuals, " + expsize + " expected size " + actsize + " actual size. NOT OK");
			} else {
				System.out.println(indir + "GenotypeMatrix.dat: " + nrsnps + " snps, " + nrinds + " individuals, " + expsize + " expected size " + actsize + " actual size. OK");
			}
		}
		
		if (Gpio.exists(indir + "ImputedDosageMatrix.dat")) {
			long expsize = ((long) nrinds * 1) * (long) nrsnps;
			long actsize = Gpio.getFileSize(indir + "ImputedDosageMatrix.dat");
			if (expsize != actsize) {
				System.out.println(indir + "ImputedDosageMatrix.dat: " + nrsnps + " snps, " + nrinds + " individuals, " + expsize + " expected size " + actsize + " actual size. NOT OK");
			} else {
				System.out.println(indir + "ImputedDosageMatrix.dat: " + nrsnps + " snps, " + nrinds + " individuals, " + expsize + " expected size " + actsize + " actual size. OK");
			}
		}
		
	}
}
