package nl.harmjanwestra.playground.conv;

import nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;

public class TriTyperDatasetFilter {

    public static void main(String[] args) {
        TriTyperDatasetFilter d = new TriTyperDatasetFilter();

        if (args.length < 2) {
            System.out.println("Usage: indir outdir [maf] [cr] [hwep]");
        } else {
            try {

                double maf = 0.01;
                double hwep = 0.0001;
                double cr = 0.95;

                if (args.length > 2) {
                    maf = Double.parseDouble(args[2]);
                }
                if (args.length > 3) {
                    cr = Double.parseDouble(args[3]);
                }
                if (args.length > 4) {
                    hwep = Double.parseDouble(args[4]);
                }

                d.run(args[0], args[1], maf, cr, hwep);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public void run(String indir, String outdir, double maf, double cr, double hwep) throws IOException {

        TriTyperGenotypeData d = new TriTyperGenotypeData(indir);
        SNPLoader loader = d.createSNPLoader(10000);
        String[] snps = d.getSNPs();
        String[] inds = d.getIndividuals();
        ArrayList<Integer> snpsToKeep = new ArrayList<Integer>();

        ProgressBar pb = new ProgressBar(snps.length, "Checking SNPs.");
        for (int s = 0; s < snps.length; s++) {
            SNP obj = d.getSNPObject(s);
            loader.loadGenotypes(obj);
            if (obj.getMAF() > maf && obj.getCR() > cr && obj.getHWEP() > hwep) {
                snpsToKeep.add(s);
            }
            obj.clearGenotypes();
            obj = null;
            pb.iterate();
        }
        pb.close();
        System.out.println(snpsToKeep.size() + " snps to keep");
        loader.close();

        BinaryFile bf = new BinaryFile(outdir + "GenotypeMatrix.dat", BinaryFile.W);
        Gpio.copyFile(indir + "Individuals.txt", outdir + "Individuals.txt");
        Gpio.copyFile(indir + "PhenotypeInformation.txt", outdir + "PhenotypeInformation.txt");

        RandomAccessFile rf = new RandomAccessFile(indir + "GenotypeMatrix.dat", "r");
        TextFile snpout = new TextFile(outdir + "SNPs.txt.gz", TextFile.W);
        TextFile snpmapout = new TextFile(outdir + "SNPMappings.txt.gz", TextFile.W);
        int nrBytesToRead = inds.length * 2;
        pb = new ProgressBar(snpsToKeep.size(), "Writing SNPs.");
        for (int i = 0; i < snpsToKeep.size(); i++) {
            long seekLoc = (long) i * (long) nrBytesToRead;
            long seekEnd = seekLoc + (long) (inds.length * 2);

            rf.seek(seekLoc);
            byte[] data = new byte[(int) (seekEnd - seekLoc)];
            rf.read(data);
            bf.write(data);

            snpout.writeln(snps[i]);
            snpmapout.writeln(d.getChr(i) + "\t" + d.getChrPos(i) + "\t" + snps[i]);
            pb.iterate();
        }
        pb.close();
        rf.close();
        bf.close();
        snpout.close();
        snpmapout.close();

    }
}
