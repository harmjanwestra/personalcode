package nl.harmjanwestra.playground.conv;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.bin.RandomAccessFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class TriTyperDatasetSorter {


    public static void main(String[] args) {


        TriTyperDatasetSorter sorter = new TriTyperDatasetSorter();
        if (args.length < 4) {
            System.out.println("Usage: inputfolder snps.txt.gz snpmappings.txt.gz outdir [dropunmapping:false] [dropduplicates:false]");
        } else {
            String folder = args[0];
            String snps = args[1];
            String snpmap = args[2];
            String outdir = args[3];
            boolean dropunmapping = false;
            if (args.length > 4) {
                dropunmapping = Boolean.parseBoolean(args[4]);
            }
            boolean dropduplicates = false;
            if (args.length > 5) {
                dropduplicates = Boolean.parseBoolean(args[5]);
            }
            try {
                sorter.run(folder, snps, snpmap, dropunmapping, dropduplicates, outdir);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }


    }


    private class SNPObj {
        String id;
        long gtbytepos = -1;
        long impbytepos = -1;
        short chr;
        int pos;

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            SNPObj snpObj = (SNPObj) o;
            return chr == snpObj.chr &&
                    pos == snpObj.pos;
        }

        @Override
        public int hashCode() {
            return Objects.hash(chr, pos);
        }
    }

    public void run(String inputDir, String snps, String snpmap, boolean dropunmapping, boolean dropduplicates, String outputDir) throws IOException {

        System.out.println("Input: " + inputDir);
        System.out.println("SNPs: " + snps);
        System.out.println("SNPMappings: " + snpmap);
        System.out.println("Outdir: " + outputDir);
        System.out.println("Dropunmapping: " + dropunmapping);
        System.out.println("Dropduplicates: " + dropduplicates);


        Gpio.createDir(outputDir);


        // read SNPset
        HashMap<String, SNPObj> snpToObj = new HashMap<>();
        TextFile tf = new TextFile(snpmap, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        System.out.println("Reading " + snpmap);
        int lctr = 0;
        while (elems != null) {

            int chr = ChrAnnotation.parseChr(elems[0]);
            int pos = Integer.parseInt(elems[1]);
            String id = new String(elems[2]);

            if (!dropunmapping || (chr > 0 && chr < 26)) {
                SNPObj snpObj = new SNPObj();
                snpObj.chr = (short) chr;
                snpObj.pos = pos;
                snpObj.id = id;
                snpToObj.put(id, snpObj);
            }

            lctr++;
            if (lctr % 1000000 == 0) {
                System.out.print(lctr + " lines read\r");
            }

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println();
        System.out.println(snpToObj.size() + " ids in file");

        tf = new TextFile(inputDir + "/Individuals.txt", TextFile.R);
        String ind = tf.readLine();
        ArrayList<String> individuals = new ArrayList<>();
        while (ind != null) {
            if (ind.length() > 0) {
                individuals.add(ind);
            }
            ind = tf.readLine();
        }
        tf.close();
        System.out.println(individuals.size() + " individuals found");
        tf = new TextFile(snps, TextFile.R);
        String ln = tf.readLine();
        long curpos = 0;
        long curimppos = 0;
        System.out.println("Reading: " + snps);
        int nrinds = individuals.size();
        int found = 0;
        int snpctr = 0;
        while (ln != null) {

            String id = ln.trim();
            SNPObj obj = snpToObj.get(id);
            if (obj != null) {
                obj.gtbytepos = curpos;
                obj.impbytepos = curimppos;
                found++;
            }

            curpos = curpos + (nrinds * 2);
            curimppos = curimppos + nrinds;
            snpctr++;
            if (snpctr % 1000000 == 0) {
                System.out.print(snpctr + " lines read. " + found + " SNPs found\r");
            }
            ln = tf.readLine();
        }
        tf.close();
        System.out.println();
        System.out.println(found + " snps found out of " + snpctr);

        // check if matrices have correct length
        long explen = ((long) individuals.size() * 2) * (long) snpctr;
        File f = new File(inputDir + "GenotypeMatrix.dat");
        long actuallen = f.length();
        boolean ok = true;
        if (actuallen != explen) {
            System.err.println("Error: expected " + explen + " but found " + actuallen + " for " + f.getAbsolutePath());
            ok = false;
        }

        if (Gpio.exists(inputDir + "ImputedDosageMatrix.dat")) {
            long explendosage = ((long) individuals.size()) * (long) snpctr;
            File f2 = new File(inputDir + "ImputedDosageMatrix.dat");
            long actuallendosage = f2.length();
            if (explendosage != actuallendosage) {
                System.err.println("Error: expected " + explendosage + " but found " + actuallendosage + " for " + f2.getAbsolutePath());
                ok = false;
            }
        }

        if (!ok) {
            System.exit(-1);
        }

        ArrayList<SNPObj> objs = new ArrayList<SNPObj>();
        for (String id : snpToObj.keySet()) {
            SNPObj obj = snpToObj.get(id);
            if (obj.impbytepos > -1 && obj.gtbytepos > -1) {
                objs.add(obj);
            }
        }
        System.out.println(objs.size() + " SNPs after filter");
        SNPObj[] objarr = objs.toArray(new SNPObj[0]);
        Arrays.parallelSort(objarr, new SNPObjSorter());


        System.out.println("writing to " + outputDir);
        Gpio.copyFile(inputDir + "PhenotypeInformation.txt", outputDir + "PhenotypeInformation.txt");
        Gpio.copyFile(inputDir + "Individuals.txt", outputDir + "Individuals.txt");
        TextFile snpout = new TextFile(outputDir + "SNPs.txt.gz", TextFile.W);
        TextFile snpmapout = new TextFile(outputDir + "SNPMappings.txt.gz", TextFile.W);
        for (int ctr = 0; ctr < objarr.length; ctr++) {
            SNPObj obj = objarr[ctr];
            snpout.writeln(obj.id);
            snpmapout.writeln(obj.chr + "\t" + obj.pos + "\t" + obj.id);
        }
        snpout.close();
        snpmapout.close();
        System.out.println("SNPMap and SNPs.txt.gz written");


        int buffersize = 10000;
        System.out.println("Reading genotypes");
        BinaryFile bf = new BinaryFile(outputDir + "GenotypeMatrix.dat", BinaryFile.W);
        RandomAccessFile raf = new RandomAccessFile(inputDir + "GenotypeMatrix.dat", RandomAccessFile.READ);
        byte[][] buffer = new byte[buffersize][nrinds * 2];
        int bctr = 0;
        ProgressBar pb = new ProgressBar(objarr.length, "Reading SNPs");
        for (int ctr = 0; ctr < objarr.length; ctr++) {
            SNPObj obj = objarr[ctr];
            raf.seek(obj.gtbytepos);
            raf.read(buffer[bctr]);
            bctr++;
            if (bctr == buffersize) {
                for (int c = 0; c < bctr; c++) {
                    bf.write(buffer[c]);
                }
                bctr = 0;
            }
            pb.iterate();
        }

        // write remaining buffer
        for (int c = 0; c < bctr; c++) {
            bf.write(buffer[c]);
        }
        pb.close();
        raf.close();
        bf.close();

        if (Gpio.exists(inputDir + "ImputedDosageMatrix.dat")) {
            buffer = new byte[buffersize][nrinds];
            BinaryFile bfdosage = new BinaryFile(outputDir + "ImputedDosageMatrix.dat", BinaryFile.W);
            RandomAccessFile rafd = new RandomAccessFile(inputDir + "ImputedDosageMatrix.dat", RandomAccessFile.READ);
            bctr = 0;
            pb = new ProgressBar(objarr.length, "Reading SNP dosagess");
            pb.set(0);
            for (int ctr = 0; ctr < objarr.length; ctr++) {
                SNPObj obj = objarr[ctr];
                rafd.seek(obj.impbytepos);
                rafd.read(buffer[bctr]);
                bctr++;
                if (bctr == buffersize) {
                    for (int c = 0; c < bctr; c++) {
                        bfdosage.write(buffer[c]);
                    }
                    bctr = 0;
                }
                pb.iterate();
            }

            // write remaining bufferb
            for (int c = 0; c < bctr; c++) {
                bfdosage.write(buffer[c]);
            }
            pb.close();

            bfdosage.close();
            rafd.close();
        }

        System.out.println("Done. Have a nice day.");

    }

    private class SNPObjSorter implements Comparator<SNPObj> {
        @Override
        public int compare(SNPObj snpObj, SNPObj t1) {
            if (snpObj.equals(t1)) {
                return 0;
            } else {
                if (snpObj.chr == t1.chr) {
                    return Integer.compare(snpObj.pos, t1.pos);
                } else {
                    return Integer.compare(snpObj.chr, t1.chr);
                }
            }
        }
    }
}
