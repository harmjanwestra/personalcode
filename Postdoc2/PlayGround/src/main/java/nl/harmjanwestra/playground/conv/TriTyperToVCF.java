package nl.harmjanwestra.playground.conv;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class TriTyperToVCF {

    public static void main(String[] args) {
        TriTyperToVCF d = new TriTyperToVCF();

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
        System.out.println("TriTyper VCF converter");
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

        TextFile vcfout = new TextFile(out + "ConvertedFromTriTyper.vcf.gz", TextFile.W);
        vcfout.writeln("##fileformat=VCFv4.2");

        String header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + Strings.concat(individuals, Strings.tab);
        vcfout.writeln(header);

        int nrinds = individuals.length;

        byte[] a1 = new byte[nrinds];
        byte[] a2 = new byte[nrinds];

        Boolean[] inc = new Boolean[nrinds];
        Boolean[] female = new Boolean[nrinds];
        for (int i = 0; i < nrinds; i++) {
            inc[i] = true;
        }


        HashMap<String, int[]> mappings = new HashMap<String, int[]>();
        TextFile tfq = new TextFile(loc + "SNPMappings.txt.gz", TextFile.R);
        String[] elems = tfq.readLineElems(TextFile.tab);
        System.out.println("Reading SNP map: " + loc + "SNPMappings.txt.gz");
        int lnctr = 0;
        while (elems != null) {
            int chr = Chromosome.parseChr(elems[0]).getNumber();
            int pos = Integer.parseInt(elems[1]);
            String rs = elems[2];

            int[] chrpos = new int[]{chr, pos};
            mappings.put(rs, chrpos);
            lnctr++;
            if (lnctr % 10000 == 0) {
                System.out.print(lnctr + " read sofar.\r");
            }
            elems = tfq.readLineElems(TextFile.tab);
        }
        tfq.close();
        System.out.println();
        System.out.println(mappings.size() + " mappings read.");

        umcg.genetica.io.text.TextFile tf2 = new umcg.genetica.io.text.TextFile(loc + "SNPs.txt.gz", umcg.genetica.io.text.TextFile.R);
        String ln = tf2.readLine();
        int ctr = 0;
        int written = 0;
        BinaryFile bf = new BinaryFile(loc + "GenotypeMatrix.dat", BinaryFile.R, 32 * 1024);
        TextFile log = new TextFile(out + "ConversionLog.txt.gz", TextFile.W);
        log.writeln("SNP\tNrCalled\tCallRate\tMAF\tHWE-P");
        while (ln != null) {
            bf.read(a1);
            bf.read(a2);
            SNP snp = new SNP();
            snp.setAlleles(a1, a2, inc, female);
//            System.out.println(ln + "\t" + snp.getMAF() + "\t" + snp.nrCalled + "\t" + snp.getCR() + "\t" + snp.getHWEP());
            log.writeln(ln + "\t" + snp.getNrCalled() + "\t" + snp.getCR() + "\t" + snp.getMAF() + "\t" + snp.getHWEP());
            if (snp.getMAF() > maf && snp.getCR() > cr && snp.getHWEP() > hwep) {
                // write

                byte al1 = snp.getAlleles()[0];
                byte al2 = snp.getAlleles()[1];
                String[] outputgt = new String[nrinds];
                byte[] gt = snp.getGenotypes();
                for (int i = 0; i < nrinds; i++) {
                    byte indgt = gt[i];
                    if (indgt == -1) {
                        outputgt[i] = "./.";
                    } else if (indgt == 0) {
                        outputgt[i] = "0/0";
                    } else if (indgt == 1) {
                        outputgt[i] = "0/1";
                    } else if (indgt == 2) {
                        outputgt[i] = "1/1";
                    }
                }

                int[] chrpos = mappings.get(ln);
                if (chrpos != null) {
                    String outln = chrpos[0] + "\t" + chrpos[1] + "\t" + ln + "\t" + BaseAnnot.toString(al1) + "\t" + BaseAnnot.toString(al2) + "\t.\t.\t.\tGT\t" + Strings.concat(outputgt, Strings.tab);
                    vcfout.writeln(outln);
                } else {
                    String outln = 0 + "\t" + 0 + "\t" + ln + "\t" + BaseAnnot.toString(al1) + "\t" + BaseAnnot.toString(al2) + "\t.\t.\t.\tGT\t" + Strings.concat(outputgt, Strings.tab);
                    vcfout.writeln(outln);
                }

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
        log.close();
        vcfout.close();

        System.out.println();
        System.out.print(ctr + " lines parsed total. " + written + " snps written");

    }
}
