package nl.harmjanwestra.playground.conv;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.trityper.SNP;

import java.io.IOException;

public class TriTyperDatasetMAFFilter {

    public static void main(String[] args) {
        TriTyperDatasetMAFFilter d = new TriTyperDatasetMAFFilter();
//        args = new String[]{
//                "D:\\Work\\ampad\\dl\\TriTyper-merged\\",
//                "D:\\Work\\ampad\\dl\\TriTyper-merged\\filter\\"
//        };
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

    public void run(String loc, String out, double maf, double cr, double hwep) throws IOException {
        System.out.println("TriTyper MAF filter");
        System.out.println("In: " + loc);
        System.out.println("Out: " + out);
        System.out.println("MAF: " + maf);
        System.out.println("CR: " + cr);
        System.out.println("HWEP: " + hwep);
        umcg.genetica.io.text.TextFile tf = new umcg.genetica.io.text.TextFile(loc + "Individuals.txt", umcg.genetica.io.text.TextFile.R);
        String[] individuals = tf.readAsArray();
        tf.close();
        System.out.println(individuals.length + " individuals");
        int nrBytesToRead = individuals.length * 2;

        BinaryFile bf = new BinaryFile(loc + "GenotypeMatrix.dat", BinaryFile.R, 32 * 1024);
        BinaryFile bfo = new BinaryFile(out + "GenotypeMatrix.dat", BinaryFile.W, 32 * 1024);
        int nrinds = individuals.length;

        byte[] a1 = new byte[nrinds];
        byte[] a2 = new byte[nrinds];

        Boolean[] inc = new Boolean[nrinds];
        Boolean[] female = new Boolean[nrinds];
        for (int i = 0; i < nrinds; i++) {
            inc[i] = true;
        }


        umcg.genetica.io.text.TextFile tf2 = new umcg.genetica.io.text.TextFile(loc + "SNPs.txt.gz", umcg.genetica.io.text.TextFile.R);
        umcg.genetica.io.text.TextFile tfsnpout = new umcg.genetica.io.text.TextFile(out + "SNPs.txt.gz", umcg.genetica.io.text.TextFile.W);
        String ln = tf2.readLine();
        int ctr = 0;
        int written = 0;
        while (ln != null) {
            bf.read(a1);
            bf.read(a2);
            SNP snp = new SNP();
            snp.setAlleles(a1, a2, inc, female);
            // System.out.println(ln + "\t" + snp.getMAF() + "\t" + snp.nrCalled + "\t" + snp.getCR() + "\t" + snp.getHWEP());
            if (snp.getMAF() > maf && snp.getCR() > cr && snp.getHWEP() > hwep) {
                // write
                tfsnpout.writeln(ln);
                bfo.write(a1);
                bfo.write(a2);
                written++;
            }
            ctr++;
            if (ctr % 10000 == 0) {
                System.out.print("\r" + ctr + " lines parsed sofar. " + written + " snps written");
            }
            ln = tf2.readLine();
        }

        tf2.close();
        bf.close();
        bfo.close();
        tfsnpout.close();

        Gpio.copyFile(loc + "Individuals.txt", out + "Individuals.txt");
        Gpio.copyFile(loc + "PhenotypeInformation.txt", out + "PhenotypeInformation.txt");
        Gpio.copyFile(loc + "SNPMappings.txt.gz", out + "SNPMappings.txt.gz");
        System.out.println();
        System.out.print(ctr + " lines parsed total. " + written + " snps written");

    }
}
