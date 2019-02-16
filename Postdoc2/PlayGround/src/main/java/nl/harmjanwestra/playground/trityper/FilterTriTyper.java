package nl.harmjanwestra.playground.trityper;


import nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

import java.io.IOException;
import java.io.RandomAccessFile;

public class FilterTriTyper {
	
	public static void main(String[] args) {
		String loc = "D:\\Sync\\SyncThing\\geuvadis\\genotypes-EUR\\";
		String out = "D:\\Sync\\SyncThing\\geuvadis\\genotypes-EUR-10perc\\";
		double perc = 0.1;
		
		FilterTriTyper t = new FilterTriTyper();
		try {
			t.filter(loc, out, perc);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void filter(String loc, String out, double perc) throws IOException {
		
		TriTyperGenotypeData ds = new TriTyperGenotypeData(loc);
		RandomAccessFile bf = new RandomAccessFile(loc + "GenotypeMatrix.dat", "r");
		
		BinaryFile f = new BinaryFile(out + "GenotypeMatrix.dat", BinaryFile.W);
		int nrBytesToRead = ds.getIndividuals().length * 2;
		String[] snps = ds.getSNPs();
		byte[] data = new byte[nrBytesToRead];
		TextFile tf = new TextFile(out + "SNPs.txt.gz", TextFile.W);
		TextFile tf2 = new TextFile(out + "SNPMappings.txt.gz", TextFile.W);
		ProgressBar pb = new ProgressBar(snps.length);
		for (int i = 0; i < snps.length; i++) {
			if (Math.random() < perc) {
				SNP snp = ds.getSNPObject(i);
				long seekLoc = (long) snp.getId() * (long) nrBytesToRead;
				bf.seek(seekLoc);
				bf.read(data);
				f.write(data);
				tf.writeln(snps[i]);
				tf2.writeln(snp.getChr() + "\t" + snp.getChrPos() + "\t" + snp.getName());
			}
			pb.set(i);
		}
		pb.close();
		tf.close();
		tf2.close();
		f.close();
		bf.close();
		
		Gpio.copyFile(loc + "Individuals.txt", out + "Individuals.txt");
		Gpio.copyFile(loc + "PhenotypeInformation.txt", out + "PhenotypeInformation.txt");
		
		
	}
	
	
}
