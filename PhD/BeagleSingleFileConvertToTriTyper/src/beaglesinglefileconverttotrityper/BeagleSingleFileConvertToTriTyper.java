/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package beaglesinglefileconverttotrityper;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.WGAFileMatrixImputedDosage;

/**
 *
 * @author harmjan
 */
public class BeagleSingleFileConvertToTriTyper {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        if (args.length < 2) {
            System.out.println("Usage: infile outdir");
        } else {
            try {
                BeagleSingleFileConvertToTriTyper.run(args[0], args[1]);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public static void run(String file, String outputDir) throws IOException {
        if(!outputDir.endsWith("/")){
            outputDir += "/";
        }
        Gpio.createDir(outputDir);
        
        TextFile in = new TextFile(file, TextFile.R);
        int prevLine = -1;
        int line = 0;
        int numValues = 0;


        HashMap<String, Integer> hashSNP = new HashMap<String, Integer>();
        ArrayList<String> ArrayListSNP = new ArrayList<String>();
        ArrayList<String> ArrayListSNPMappings = new ArrayList<String>();
        ArrayList<String> ArrayListInd = new ArrayList<String>();
        HashMap<String, Integer> hashInd = new HashMap<String, Integer>();
        long nrSNPsAvailable = 0;
        String[] data;


        String str = in.readLine(); // header
        data = str.split(" ");

        for (int c = 3; c < data.length; c++) {
            if (!hashInd.containsKey(data[c])) {
                hashInd.put(data[c], ArrayListInd.size());
                ArrayListInd.add(data[c]);
                System.out.println("Found new individual:" + data[c]);

                // System.out.println(data[c]);
            }
        }
        
        str = in.readLine(); // first line
        while (str != null) {
            while (str.contains("  ")) {
                str = str.replace("  ", " ");
            }
            data = str.split(" ");

            numValues = (data.length - 3) / 3;

            String snp = new String(data[0].getBytes());
            if (!hashSNP.containsKey(snp)) {
                String snpPos = "1";
                String snpMapping = 0 + "\t" + snpPos + "\t" + snp;
                ArrayListSNPMappings.add(snpMapping);
                hashSNP.put(snp, ArrayListSNP.size());
                ArrayListSNP.add(snp);
                nrSNPsAvailable++;

            }

            if (line % 10000 == 0 && line > prevLine) {
                System.out.print(".");
                prevLine = line;
            }
            line++;
            str = in.readLine();
        }
        System.out.println("");
        System.out.println("Number of SNPs parsed so far:\t" + nrSNPsAvailable + " for " + numValues + " samples");

        System.out.println("");
        in.close();


        System.out.println("Number of individuals parsed:\t" + ArrayListInd.size());
        System.out.println("\nWriting SNP mappings to file:");

        BufferedWriter snpmappingsout = new BufferedWriter(new FileWriter(outputDir + "SNPMappings.txt"));
        for (int snp = 0; snp < ArrayListSNPMappings.size(); snp++) {
            snpmappingsout.write(((String) ArrayListSNPMappings.get(snp)) + "\n");
            snpmappingsout.flush();
            if (snp % 2000 == 1999) {
                System.out.print(".");
            }
        }
        System.out.println("");
        snpmappingsout.close();


        System.out.println("\nWriting marker definition to file:");

        BufferedWriter snpstxt = new BufferedWriter(new FileWriter(outputDir + "SNPs.txt"));
        for (int snp = 0; snp < ArrayListSNP.size(); snp++) {
            snpstxt.write(((String) ArrayListSNP.get(snp)) + "\n");
            snpstxt.flush();
            if (snp % 2000 == 1999) {
                System.out.print(".");
            }
        }
        System.out.println("");
        snpstxt.close();


        System.out.println("\nWriting individuals to file:");

        BufferedWriter outInd = new BufferedWriter(new FileWriter(outputDir + "Individuals.txt"));

        BufferedWriter outPhe = new BufferedWriter(new FileWriter(outputDir + "PhenotypeInformation.txt"));

        for (int ind = 0; ind < ArrayListInd.size(); ind++) {
            outInd.write(((String) ArrayListInd.get(ind)) + "\n");

            outPhe.write(((String) ArrayListInd.get(ind)) + "\tcontrol\tinclude\tunknown" + "\n");

            outInd.flush();
            if (ind % 5 == 4) {
                System.out.print(".");
            }
        }

        System.out.println("");
        outInd.close();
        outPhe.close();

        int nrSNPs = (int) nrSNPsAvailable;
        int nrSamples = ArrayListInd.size();
        WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(nrSNPs, nrSamples, new File(outputDir + "GenotypeMatrix.dat"), false);
        WGAFileMatrixImputedDosage matrixImputedDosage = new WGAFileMatrixImputedDosage(nrSNPs, nrSamples, new File(outputDir + "/ImputedDosageMatrix.dat"), false);




        System.out.print("Processing file:\t" + file);

        // BufferedReader in = new BufferedReader(new FileReader(new File(fileName)));
        //BufferedReader in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fileName))));
        in.open();

        str = in.readLine();

        prevLine = -1;
        line = 0;
        str = in.readLine();
        while (str != null) {
            while (str.contains("  ")) {
                str = str.replace("  ", " ");
            }
            data = str.split(" ");
            String snp = new String(data[0].getBytes());
            int snpIndex = ((Integer) hashSNP.get(snp)).intValue();
            byte[] allele1 = new byte[nrSamples];
            byte[] allele2 = new byte[nrSamples];
            byte[] alleles = new byte[2];
            alleles[0] = data[1].getBytes()[0];
            alleles[1] = data[2].getBytes()[0];
            byte[] dosage = new byte[nrSamples];
            for (int sample = 0; sample < nrSamples; sample++) {
                // AB BB
                double dosageValue = Double.parseDouble(data[sample * 3 + 4]) * 1d + Double.parseDouble(data[sample * 3 + 5]) * 2d;

                int dosageInt = (int) Math.round(dosageValue * 100d);
                byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                if (dosageInt < 0 || dosageInt > 200) {
                    System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpIndex + "\t" + data[sample * 3 + 3] + "-" + data[sample * 3 + 4] + "-" + data[sample * 3 + 5]);
                } else {
                    dosage[sample] = (byte) dosageByte;
                }
                if (dosageValue < 0.5) {
                    allele1[sample] = alleles[0];
                    allele2[sample] = alleles[0];
                } else {
                    if (dosageValue > 1.5) {
                        allele1[sample] = alleles[1];
                        allele2[sample] = alleles[1];
                    } else {
                        allele1[sample] = alleles[0];
                        allele2[sample] = alleles[1];
                    }
                }
                fileMatrixGenotype.setAllele1(snpIndex, sample, allele1);
                fileMatrixGenotype.setAllele2(snpIndex, sample, allele2);
                matrixImputedDosage.setDosage(snpIndex, sample, dosage);
            }

            if (line % 10000 == 0 && line > prevLine) {
                System.out.print(".");
                prevLine = line;
            }
            line++;
            str = in.readLine();
        }
        in.close();

        System.out.println("");


        fileMatrixGenotype.close();

        matrixImputedDosage.close();
    }
}
