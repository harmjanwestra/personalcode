package nl.harmjanwestra.playground.conv;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.playground.legacy.vcf.VCFGenotypeData;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

public class ConvertVCFToTT {


    public static void main(String[] args) {


        ConvertVCFToTT v = new ConvertVCFToTT();
//		String[] allowedVQSROrFilter = new String[]{
//				"PASS", "GENOTYPED",
//				"GENOTYPED_ONLY",
//				"VQSRTrancheINDEL99.90to99.95",
//				"VQSRTrancheINDEL99.00to100.00+",
//				"VQSRTrancheINDEL99.00to100.00",
//				"VQSRTrancheINDEL99.95to100.00+",
//				"VQSRTrancheINDEL99.95to100.00",
//				"VQSRTrancheSNP99.80to100.00",
//				"VQSRTrancheSNP99.80to99.90",
//				"VQSRTrancheSNP99.90to99.95",
//				"VQSRTrancheSNP99.95to100.00+",
//				"VQSRTrancheSNP99.95to100.00",
//				"VQSRTrancheSNP99.80to100.00",
//				"VQSRTrancheSNP99.80to100.00+",
//				"VQSRTrancheSNP99.80to100.00"
//		};
//		HashSet<String> allowedFilter = new HashSet<String>();
//		allowedFilter.addAll(Arrays.asList(allowedVQSROrFilter));

        boolean writegenotypes = true;
        boolean parsevcf = true;
        boolean useAD = true;

        if (args.length < 2) {
            System.out.println("Usage: inputvcf outputfolder [maf:0.001] [missingness:0.05] [minAD:10] [ablower:0.2]" +
                    " [abupper:0.8] [gq:20] [impqual:0.3] [inbreedingcoeff:-0.3] [includeIndels:true] [applyInbreedingFilterToSNPs:true] [allowedInfoOrVQSRStringsSemicolonSep]");
        } else {

//            try {
//                v.ConvertWithoutFilter(args[0], args[1], 0.01, 0.5, allowedFilter);
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//
//            System.exit(-1);
            try {
                double maf = 0.001;
                if (args.length > 2) {
                    maf = Double.parseDouble(args[2]);
                }


                double missingness = 0.15;
                if (args.length > 3) {
                    missingness = Double.parseDouble(args[3]);
                }

                Integer minad = 10;
                if (args.length > 4) {
                    minad = Integer.parseInt(args[4]);
                }
                double ABlower = 0.2;
                if (args.length > 5) {
                    ABlower = Double.parseDouble(args[5]);
                }
                double ABupper = 0.8;
                if (args.length > 6) {
                    ABupper = Double.parseDouble(args[6]);
                }
                double GQ = 20;
                if (args.length > 7) {
                    GQ = Double.parseDouble(args[7]);
                }
                double impqualthreshold = 0.3;
                if (args.length > 8) {
                    impqualthreshold = Double.parseDouble(args[8]);
                }
                double inbreedingCoeff = -0.3;
                if (args.length > 9) {
                    inbreedingCoeff = Double.parseDouble(args[9]);
                }

                boolean allowIndels = true;
                if (args.length > 10) {
                    allowIndels = Boolean.parseBoolean(args[10]);
                }

                boolean applyInbreedingFilterToSNPs = true;
                if (args.length > 11) {
                    applyInbreedingFilterToSNPs = Boolean.parseBoolean(args[10]);
                }

                HashSet<String> allowedFilter = null;
                if (args.length > 12) {

                    String[] filterelems = args[12].split(";");
                    allowedFilter = new HashSet<>();
                    for (String s : filterelems) {
                        allowedFilter.add(s);
                        System.out.println("Allowing variants with filter str: " + s);
                    }
                }


                v.run(args[0],
                        args[1],
                        parsevcf,
                        writegenotypes,
                        maf, missingness,
                        minad,
                        ABlower,
                        ABupper,
                        GQ,
                        allowedFilter,
                        inbreedingCoeff,
                        impqualthreshold,
                        allowIndels,
                        applyInbreedingFilterToSNPs);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }


    }

    public void run(String vcffile, String output, boolean parseVCF, boolean writegenotypes,
                    double mafthreshold, double missingnessthreshold, Integer minAD,
                    double ABlower, double ABupper, double GQ,
                    HashSet<String> allowedVQSROrINFO, double inbreedingCoeff,
                    double impqualthreshold, boolean allowIndels, boolean applyInbreedingFilterToSNPs) throws IOException {

        System.out.println("in: " + vcffile);
        System.out.println("out: " + output);
        System.out.println("maf: " + mafthreshold);
        System.out.println("missingness: " + missingnessthreshold);
        System.out.println("minAD: " + minAD);
        System.out.println("ablower: " + ABlower);
        System.out.println("abupper: " + ABupper);
        System.out.println("GQ: " + GQ);
        System.out.println("Inbreeding: " + inbreedingCoeff);
        System.out.println("ImpQual: " + impqualthreshold);
        System.out.println("AllowIndels: " + allowIndels);
        System.out.println("ApplyInbreedingCoeffToSNPs: " + applyInbreedingFilterToSNPs);

        String folder = output;
        Gpio.createDir(folder);
        File genotypeDataFile1 = new File(folder, "GenotypeMatrix.dat");
        File imputedDosageDataFile1 = new File(folder, "ImputedDosageMatrix.dat");

        File snpFile = new File(folder, "SNPs.txt.gz");

        File snpFile2 = null;
        File genotypeDataFile2 = null;
        File imputedDosageDataFile2 = null;
        if (allowIndels) {
            genotypeDataFile2 = new File(folder, "GenotypeMatrixV2.dat");
            imputedDosageDataFile2 = new File(folder, "ImputedDosageMatrixV2.dat");
            snpFile2 = new File(folder, "SNPsV2.txt.gz");
        }

        File snpMapFile = new File(folder, "SNPMappings.txt.gz");
        File individualFile = new File(folder, "Individuals.txt");
        File phenotypeAnnotationFile = new File(folder, "PhenotypeInformation.txt");
        File logfile = new File(folder, "ConversionLog.txt.gz");


        if (!parseVCF) {
            System.exit(0);
        }

        BufferedOutputStream genotypeDataFileWriter1 = null;
        BufferedOutputStream genotypeDosageDataFileWriter1 = null;

        BufferedOutputStream genotypeDataFileWriter2 = null;
        BufferedOutputStream genotypeDosageDataFileWriter2 = null;

        if (writegenotypes) {
            genotypeDataFileWriter1 = new BufferedOutputStream(new FileOutputStream(genotypeDataFile1), 32 * 1024);
            genotypeDosageDataFileWriter1 = new BufferedOutputStream(new FileOutputStream(imputedDosageDataFile1), 32 * 1024);
            if (allowIndels) {
                genotypeDataFileWriter2 = new BufferedOutputStream(new FileOutputStream(genotypeDataFile2), 32 * 1024);
                genotypeDosageDataFileWriter2 = new BufferedOutputStream(new FileOutputStream(imputedDosageDataFile2), 32 * 1024);
            }
        }
        TextFile snpFileWriter1 = new TextFile(snpFile, TextFile.W);
        TextFile snpFileWriter2 = null;
        if (allowIndels) {
            snpFileWriter2 = new TextFile(snpFile2, TextFile.W);
        }
        TextFile snpMapFileWriter = new TextFile(snpMapFile, TextFile.W);


        String[] files;
        ArrayList<String> filestmp = new ArrayList<>();
        if (vcffile.contains("CHR")) {
            for (int i = 1; i < 25; i++) {
                String file = vcffile.replace("CHR", "" + i);

                if (i == 23) {
                    file = vcffile.replace("CHR", "X");
                } else if (i == 24) {
                    file = vcffile.replace("CHR", "Y");
                }
                if (Gpio.exists(file)) {
                    System.out.println("Found file: " + file);
                    filestmp.add(file);
                }
            }
        } else {
            filestmp.add(vcffile);
        }
        files = filestmp.toArray(new String[0]);

        int nrvarswithdosage = 0;
        int written = 0;
        int variantswithmissingness = 0;

        int nrindels = 0;
        int nrmultiallelic = 0;
        int lowmaf = 0;
        int failfilter = 0;
        int failinbreeding = 0;
        int nrlowimpqual = 0;
        ArrayList<String> samples = null;


        // check sample order
        for (int f = 0; f < files.length; f++) {
            String infile = files[f];
            System.out.println("Checking: " + infile);
            if (f == 0) {
                VCFGenotypeData d = new VCFGenotypeData(infile);
                samples = d.getSamples();
                d.close();
                BufferedWriter phenowriter = new BufferedWriter(new FileWriter(phenotypeAnnotationFile));
                TextFile indsout = new TextFile(individualFile, TextFile.W);
                for (int s = 0; s < samples.size(); s++) {
                    indsout.writeln(samples.get(s));
                    phenowriter.write(samples.get(s) + "\tunknown\tinclude\tunknown\n");
                }
                indsout.close();
                phenowriter.close();
                System.out.println(samples.size() + " samples detected.");
            } else {
                VCFGenotypeData d = new VCFGenotypeData(infile);
                ArrayList<String> sampletmp = d.getSamples();

                d.close();
                if (samples.size() == sampletmp.size()) {
                    for (int i = 0; i < samples.size(); i++) {
                        if (!samples.get(i).equals(sampletmp.get(i))) {
                            System.out.println("Issue with: " + infile);
                            System.out.println("Samples are not in the same order.");
                            System.exit(-1);
                        }
                    }
                } else {
                    System.out.println("Issue with: " + infile);
                    System.out.println("Sample sizes not equal to first file: " + sampletmp.size() + " found, " + samples.size() + " expected.");
                    System.exit(-1);
                }

            }
        }

        /*

         */


        TextFile logout = new TextFile(logfile, TextFile.W);
        String header = "var\talleles\tisIndel\tisMultiAllelic\tMAF\tPassFilter\tPassInbreeding\tHashHighMissignness\tHasLowMaf\tMissingness\tnrMissing\tAvgDepth\tVarIncluded";
        logout.writeln(header);

        int totallines = 0;
        for (int f = 0; f < files.length; f++) {

            String infile = files[f];


            TextFile tf = new TextFile(infile, TextFile.R);
            String ln = tf.readLine();
            int ctr = 0;


            HashSet<String> warningsgiven = new HashSet<>();
            while (ln != null) {
                if (!ln.startsWith("#")) {

                    VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.ALL);

                    double missingness = 0;
                    int missing = 0;

                    int nrcalled = 0;
                    boolean varincluded = false;
                    double sumdepth = 0;

                    boolean passfilter = true;
                    boolean passinbreeding = true;
                    boolean passimpqual = true;
                    boolean hashighmissingness = false;
                    boolean haslowmaf = false;

                    if (var.getAlleles().length > 2) {
                        nrmultiallelic++;
                    } else {
                        String filter = var.getFilter();
                        if (!allowIndels && var.isIndel()) {
                            passfilter = false;
                            nrindels++;
                        }

                        if (allowedVQSROrINFO != null) {
                            if (filter.contains(";")) {
                                String[] filterelems = filter.split(";");
                                for (String e : filterelems) {
                                    if (!allowedVQSROrINFO.contains(e)) {
                                        passfilter = false;
                                    }
                                }
                            } else {
                                passfilter = (filter.equals(".") || allowedVQSROrINFO.contains(filter));
                            }
                        }

                        if (!passfilter) {
                            if (!warningsgiven.contains(filter)) {
                                System.out.println("Filter option not recognized: " + filter);
                                warningsgiven.add(filter);
                            }
                            failfilter++;
                        } else if (passfilter) {

                            String inbreeding = var.getInfo().get("InbreedingCoeff");
                            Double varInbreedingCoeff = -1d;
                            if (inbreeding != null) {
                                try {
                                    varInbreedingCoeff = Double.parseDouble(inbreeding);
                                } catch (NumberFormatException e) {
                                    System.out.println("Weird Inbreeding Coeff: " + inbreeding + " for variant: " + var.toString());
                                }
                            }

                            if (applyInbreedingFilterToSNPs || var.isIndel()) {
                                passinbreeding = (inbreeding == null || varInbreedingCoeff > inbreedingCoeff);
                            }

                            if (!passinbreeding) {
                                failinbreeding++;
                            } else {

                                // check imp qual, if present
                                Double rsq = var.getImputationQualityScore();
                                if (rsq != null && rsq < impqualthreshold) {
                                    nrlowimpqual++;
                                    passimpqual = false;
                                } else {
                                    byte allele1v1 = BaseAnnot.toByte(var.getAlleles()[0]);
                                    byte allele2v1 = BaseAnnot.toByte(var.getAlleles()[1]);
                                    byte allele1v2 = 100; // BaseAnnot.toByte(var.getAlleles()[0]); // for storage of indels, ref == ascii d
                                    byte allele2v2 = 101; // BaseAnnot.toByte(var.getAlleles()[1]); // for storage of indels, alt == ascii e
                                    String ref = var.getAlleles()[0];
                                    String alt = var.getAlleles()[1];
                                    boolean isindel = var.isIndel();


                                    short[] apdepth = var.getApproximateDepth();
                                    short[][] ad = var.getAllelicDepth();

                                    DoubleMatrix2D genotypes = var.getGenotypeAllelesAsMatrix2D();
                                    DoubleMatrix2D dosages;
                                    if (var.hasPosteriorProbabilities()) {
                                        dosages = var.calculateDosageFromProbabilities(var.getPosteriorProbabilities());
                                    } else {
                                        dosages = var.getDosagesAsMatrix2D();
                                    }
//							String id = var.getId();
//							if (id.length() == 1) {
                                    String id = var.getChr() + ":" + var.getPos() + ":" + ref + "_" + alt;
//							}
                                    byte[] al11 = new byte[samples.size() * 2];
                                    byte[] al12 = new byte[samples.size() * 2];


                                    byte[] dos = new byte[samples.size()];

                                    if (var.hasImputationDosages() || var.hasPosteriorProbabilities()) {
                                        nrvarswithdosage++;
                                    }

                                    missing = 0;

                                    // check if variant has GQ
                                    short[] quals = var.getGenotypeQuals();

                                    boolean missingDP = false;
                                    if (apdepth != null) {
                                        for (int i = 0; i < samples.size(); i++) {
                                            if (apdepth[i] < 0) {
                                                missingDP = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (missingDP) {
                                        apdepth = null;
                                    }


                                    // iterate samples
                                    int nrAlt = 0;
                                    int nrTotal = 0;
                                    for (int i = 0; i < samples.size(); i++) {
                                        double a1 = genotypes.getQuick(i, 0);
                                        double a2 = genotypes.getQuick(i, 1);

                                        if (apdepth != null) {
                                            sumdepth += apdepth[i];
                                        }
                                        if (a1 == -1 || a2 == -1 || (minAD != null && apdepth != null && apdepth[i] < minAD)) {
                                            missing++;
                                        } else if (quals != null && (quals[i] > -1 && quals[i] < GQ)) {
                                            // check gq
                                            missing++;
                                        } else {
                                            if (a1 != a2) {

                                                // het.. check allelic balance if allelic depth present
                                                // ad information may not always be present..
                                                boolean ok = true;
                                                if (ad != null) {
                                                    short a1d = ad[i][0];
                                                    if (a1d > 0) { // -1 == missing
                                                        short a2d = ad[i][1];
                                                        double ab = (double) a1d / (a1d + a2d);
                                                        if (ab > ABupper || ab < ABlower) {
                                                            ok = false;
                                                        }
                                                    }
                                                }

                                                if (ok) {
                                                    // preserve phase
                                                    if (a1 == 0) {
                                                        al11[i] = allele1v1;
                                                        al11[i + samples.size()] = allele2v1;
                                                        al12[i] = allele1v2;
                                                        al12[i + samples.size()] = allele2v2;
                                                    } else {
                                                        al12[i] = allele2v2;
                                                        al12[i + samples.size()] = allele1v2;
                                                        al11[i] = allele2v1;
                                                        al11[i + samples.size()] = allele1v1;
                                                    }
                                                    nrAlt++;
                                                    nrTotal += 2;
                                                } else {
                                                    missing++;
                                                }
                                            } else {
                                                // homozygous
                                                if (a1 == 0) { // hom. reference
                                                    al11[i] = allele1v1;
                                                    al11[i + samples.size()] = allele1v1;
                                                    al12[i] = allele1v2;
                                                    al12[i + samples.size()] = allele1v2;

                                                } else { // hom. alt
                                                    al11[i] = allele2v1;
                                                    al11[i + samples.size()] = allele2v1;
                                                    al12[i] = allele2v2;
                                                    al12[i + samples.size()] = allele2v2;
                                                    nrAlt += 2;
                                                    nrTotal += 2;
                                                }
                                            }


                                            if (var.hasImputationDosages() || var.hasPosteriorProbabilities()) {
                                                double dosageval = dosages.getQuick(i, 0);
                                                dos[i] = convertDosageToByte(dosageval);
//												double v1 = var.getGenotypeProbabilies().getQuick(i, 0);
//												double v2 = var.getGenotypeProbabilies().getQuick(i, 1);
//												double v3 = var.getGenotypeProbabilies().getQuick(i, 2);
//												System.out.println(var.getChr() + ":" + var.getPos() + "\t" + a1 + "/" + a2 + "\t" + dosageval + "\t" + dos[i] + "\t" + v1 + "," + v2 + "," + v3);
                                            }
                                        }
                                    }

                                    sumdepth /= samples.size();
                                    missingness = (double) missing / samples.size();
                                    if (missingness > missingnessthreshold) {
                                        variantswithmissingness++;
                                        hashighmissingness = true;
                                    } else {

                                        // calculate MAF
                                        double maf = (double) nrAlt / nrTotal;
                                        if (maf > 0.5) {
                                            maf = 1 - maf;
                                        }
                                        if (maf < mafthreshold) {
                                            lowmaf++;
                                            haslowmaf = true;
                                        } else {

                                            if (!isindel) {
                                                snpFileWriter1.append(id);
                                                snpFileWriter1.append("\n");
                                            } else {
                                                nrindels++;
                                            }
                                            if (allowIndels) {
                                                snpFileWriter2.append(id);
                                                snpFileWriter2.append("\t");
                                                snpFileWriter2.append(Strings.concat(var.getAlleles(), Strings.comma));
                                                snpFileWriter2.append('\n');
                                            }
                                            if (!isindel) {
                                                snpMapFileWriter.append(var.getChr());
                                                snpMapFileWriter.append('\t');
                                                snpMapFileWriter.append(String.valueOf(var.getPos()));
                                                snpMapFileWriter.append('\t');
                                                snpMapFileWriter.append(id);
                                                snpMapFileWriter.append('\t');
                                                snpMapFileWriter.append(Strings.concat(var.getAlleles(), Strings.comma));
                                                snpMapFileWriter.append('\n');
                                            }
                                            if (writegenotypes) {
                                                if (!isindel) {
                                                    genotypeDataFileWriter1.write(al11);
                                                    genotypeDosageDataFileWriter1.write(dos);
                                                }
                                                if (allowIndels) {
                                                    genotypeDataFileWriter2.write(al12);
                                                    genotypeDosageDataFileWriter2.write(dos);
                                                }
                                            }

                                            written++;
                                            varincluded = true;
                                        }
                                    }
                                }

                            }
                        }
                    }

                    String incvar = "NotWritten";
                    if (varincluded) {
                        incvar = "Written";
                    }
                    String log = var.toString()
                            + "\t" + Strings.concat(var.getAlleles(), Strings.semicolon)
                            + "\t" + var.isIndel()
                            + "\t" + var.isMultiallelic()
                            + "\t" + var.getMAF()
                            + "\t" + passfilter
                            + "\t" + passinbreeding
                            + "\t" + hashighmissingness
                            + "\t" + haslowmaf
                            + "\t" + missingness
                            + "\t" + missing
                            + "\t" + sumdepth
                            + "\t" + incvar;
                    logout.writeln(log);
                }
                if (ctr % 10000 == 0) {
                    //
                    double perc = (double) written / totallines;
                    perc *= 100;
                    DecimalFormat format = new DecimalFormat("##.##");
                    String strperc = format.format(perc);
                    System.out.println(
                            "File: " + f + "\t" +
                                    ctr + "\tlines, " +
                                    totallines + "\ttotal lines, " +
                                    nrmultiallelic + " multi, " +
                                    failfilter + " filter, " +
                                    failinbreeding + " inbred, " +
                                    nrlowimpqual + " low R2, " +
                                    variantswithmissingness + " miss, " +
                                    lowmaf + " MAF, " +
                                    nrindels + " indels, " +
                                    written + " written (" + strperc + "%) - " +
                                    vcffile);
                }
                ctr++;
                totallines++;
                ln = tf.readLine();
            }
            tf.close();
        }


        logout.close();

        System.out.println();
        System.out.println("Done");
        if (writegenotypes) {
            genotypeDataFileWriter1.close();
            genotypeDosageDataFileWriter1.close();

            if (allowIndels) {
                genotypeDataFileWriter2.close();
                genotypeDosageDataFileWriter2.close();
            }
            if (nrvarswithdosage == 0) {
                System.out.println("Removing ImputedDosageMatrix.dat because there are no dosage values available");
                Gpio.delete(folder + "ImputedDosageMatrix.dat");
                Gpio.delete(folder + "ImputedDosageMatrixV2.dat");
            }
            if (nrindels == 0) {
                System.out.println("Removing GenotypeMatrixV2.dat because there are no indels");
                Gpio.delete(folder + "GenotypeMatrixV2.dat");
                Gpio.delete(folder + "SNPsV2.txt.gz");
            }
        }
        snpFileWriter1.close();
        if (allowIndels) {
            snpFileWriter2.close();
        }
        snpMapFileWriter.close();
        System.out.println(written + " written, " + nrvarswithdosage + " with dosages. " + lowmaf + " maf<" + mafthreshold + ", " + nrindels + " indels ," + nrmultiallelic + " multiallelic, " + variantswithmissingness + " with MAF>" + mafthreshold + " and high missigness " + vcffile);


    }

    private byte convertDosageToByte(double dosageValue) {
        int dosageInt = (int) Math.round(dosageValue * 100d);
        byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
        return dosageByte;
    }

    public void ConvertWithoutFilter(String vcffile, String output, double mafthreshold,
                                     double missingnessthreshold, HashSet<String> allowedVQSROrINFO) throws IOException {
        System.out.println("in: " + vcffile);
        System.out.println("out: " + output);
        System.out.println("maf: " + mafthreshold);

        String folder = output;
        File genotypeDataFile1 = new File(folder, "GenotypeMatrix.dat");
        File snpFile = new File(folder, "SNPs.txt.gz");
        File snpMapFile = new File(folder, "SNPMappings.txt.gz");
        File individualFile = new File(folder, "Individuals.txt");
        File phenotypeAnnotationFile = new File(folder, "PhenotypeInformation.txt");
        File logfile = new File(folder, "ConversionLog.txt.gz");


        BufferedOutputStream genotypeDataFileWriter1 = null;


        genotypeDataFileWriter1 = new BufferedOutputStream(new FileOutputStream(genotypeDataFile1), 32 * 1024);
        TextFile snpFileWriter1 = new TextFile(snpFile, TextFile.W);
        TextFile snpMapFileWriter = new TextFile(snpMapFile, TextFile.W);

        TextFile logout = new TextFile(logfile, TextFile.W);

        String[] files;
        ArrayList<String> filestmp = new ArrayList<>();
        if (vcffile.contains("CHR")) {
            for (int i = 1; i < 25; i++) {
                String file = vcffile.replace("CHR", "" + i);

                if (i == 23) {
                    file = vcffile.replace("CHR", "X");
                } else if (i == 24) {
                    file = vcffile.replace("CHR", "Y");
                }
                if (Gpio.exists(file)) {
                    System.out.println("Found file: " + file);
                    filestmp.add(file);
                }
            }
        } else {
            filestmp.add(vcffile);
        }
        files = filestmp.toArray(new String[0]);

        int nrvarswithdosage = 0;
        int written = 0;
        int variantswithmissingness = 0;

        int nrindels = 0;
        int nrmultiallelic = 0;
        int lowmaf = 0;
        int failfilter = 0;
        int failinbreeding = 0;
        int nrlowimpqual = 0;
        ArrayList<String> samples = null;

        // check sample order
        for (int f = 0; f < files.length; f++) {
            String infile = files[f];
            System.out.println("Checking: " + infile);
            if (f == 0) {
                VCFGenotypeData d = new VCFGenotypeData(infile);
                samples = d.getSamples();
                d.close();
                BufferedWriter phenowriter = new BufferedWriter(new FileWriter(phenotypeAnnotationFile));
                TextFile indsout = new TextFile(individualFile, TextFile.W);
                for (int s = 0; s < samples.size(); s++) {
                    indsout.writeln(samples.get(s));
                    phenowriter.write(samples.get(s) + "\tunknown\tinclude\tunknown\n");
                }
                indsout.close();
                phenowriter.close();
                System.out.println(samples.size() + " samples detected.");
            } else {
                VCFGenotypeData d = new VCFGenotypeData(infile);
                ArrayList<String> sampletmp = d.getSamples();

                d.close();
                if (samples.size() == sampletmp.size()) {
                    for (int i = 0; i < samples.size(); i++) {
                        if (!samples.get(i).equals(sampletmp.get(i))) {
                            System.out.println("Issue with: " + infile);
                            System.out.println("Samples are not in the same order.");
                            System.exit(-1);
                        }
                    }
                } else {
                    System.out.println("Issue with: " + infile);
                    System.out.println("Sample sizes not equal to first file: " + sampletmp.size() + " found, " + samples.size() + " expected.");
                    System.exit(-1);
                }

            }
        }

        for (int f = 0; f < files.length; f++) {

            String infile = files[f];


            TextFile tf = new TextFile(infile, TextFile.R);
            String ln = tf.readLine();
            int ctr = 0;

            String header = "var\talleles\tisIndel\tisMultiAllelic\tMAF\tMissingness\tnrMissing\tAvgDepth\tVarIncluded";
            logout.writeln(header);
            HashSet<String> warningsgiven = new HashSet<>();
            while (ln != null) {
                if (!ln.startsWith("#")) {

                    VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.ALL);

                    double missingness = 0;
                    int missing = 0;

                    int nrcalled = 0;
                    boolean varincluded = false;
                    double sumdepth = 0;

                    if (var.getAlleles().length > 2) {
                        nrmultiallelic++;
                    } else {
                        String filter = var.getFilter();
                        boolean passfilter = true;
                        if (filter.contains(";")) {
                            String[] filterelems = filter.split(";");
                            for (String e : filterelems) {
                                if (!allowedVQSROrINFO.contains(e)) {
                                    passfilter = false;
                                }
                            }
                        } else {
                            passfilter = (filter.equals(".") || allowedVQSROrINFO.contains(filter));
                        }

                        if (!passfilter) {
                            if (!warningsgiven.contains(filter)) {
                                System.out.println("Filter option not recognized: " + filter);
                                warningsgiven.add(filter);
                            }
                            failfilter++;
                        } else if (passfilter) {


                            byte allele1v1 = BaseAnnot.toByte(var.getAlleles()[0]);
                            byte allele2v1 = BaseAnnot.toByte(var.getAlleles()[1]);
                            byte allele1v2 = 100; // BaseAnnot.toByte(var.getAlleles()[0]);
                            byte allele2v2 = 101; // BaseAnnot.toByte(var.getAlleles()[1]);
                            String ref = var.getAlleles()[0];
                            String alt = var.getAlleles()[1];
                            boolean isindel = var.isIndel();

                            DoubleMatrix2D genotypes = var.getGenotypeAllelesAsMatrix2D();
                            DoubleMatrix2D dosages = var.getDosagesAsMatrix2D();

                            String id = var.getChr() + ":" + var.getPos() + ":" + ref + "_" + alt;
                            byte[] al11 = new byte[samples.size() * 2];
                            byte[] al12 = new byte[samples.size() * 2];

                            byte[] dos = new byte[samples.size()];

                            if (var.hasImputationDosages()) {
                                nrvarswithdosage++;
                            }

                            missing = 0;

                            // check if variant has GQ
                            boolean missingDP = false;

                            // iterate samples
                            int nrAlt = 0;
                            int nrTotal = 0;
                            for (int i = 0; i < samples.size(); i++) {
                                double a1 = genotypes.getQuick(i, 0);
                                double a2 = genotypes.getQuick(i, 1);

                                if (a1 == -1 || a2 == -1) {
                                    missing++;
                                } else {
                                    if (a1 != a2) {
                                        // het.. check allelic balance if allelic depth present
                                        // ad information may not always be present..
                                        boolean ok = true;


                                        // preserve phase
                                        if (a1 == 0) {
                                            al11[i] = allele1v1;
                                            al11[i + samples.size()] = allele2v1;
                                            al12[i] = allele1v2;
                                            al12[i + samples.size()] = allele2v2;
                                        } else {
                                            al12[i] = allele2v2;
                                            al12[i + samples.size()] = allele1v2;
                                            al11[i] = allele2v1;
                                            al11[i + samples.size()] = allele1v1;
                                        }
                                        nrAlt++;
                                        nrTotal += 2;

                                    } else {
                                        // homozygous
                                        if (a1 == 0) { // hom. reference
                                            al11[i] = allele1v1;
                                            al11[i + samples.size()] = allele1v1;
                                            al12[i] = allele1v2;
                                            al12[i + samples.size()] = allele1v2;

                                        } else { // hom. alt
                                            al11[i] = allele2v1;
                                            al11[i + samples.size()] = allele2v1;
                                            al12[i] = allele2v2;
                                            al12[i + samples.size()] = allele2v2;
                                            nrAlt += 2;
                                            nrTotal += 2;
                                        }
                                    }

                                    if (var.hasImputationDosages()) {
                                        dos[i] = convertDosageToByte(dosages.getQuick(i, 0));
                                    }
                                }
                            }

                            sumdepth /= samples.size();
                            missingness = (double) missing / samples.size();
                            if (missingness > missingnessthreshold) {
                                variantswithmissingness++;
                            } else {

                                // calculate MAF
                                double maf = (double) nrAlt / nrTotal;
                                if (maf > 0.5) {
                                    maf = 1 - maf;
                                }
                                if (maf < mafthreshold) {
                                    lowmaf++;
                                } else {

                                    if (!isindel) {
                                        snpFileWriter1.append(id);
                                        snpFileWriter1.append("\n");
                                    } else {
                                        nrindels++;
                                    }
                                    snpMapFileWriter.append(var.getChr());
                                    snpMapFileWriter.append('\t');
                                    snpMapFileWriter.append(String.valueOf(var.getPos()));
                                    snpMapFileWriter.append('\t');
                                    snpMapFileWriter.append(id);
                                    snpMapFileWriter.append('\t');
                                    snpMapFileWriter.append(Strings.concat(var.getAlleles(), Strings.comma));
                                    snpMapFileWriter.append('\n');

                                    if (!isindel) {
                                        genotypeDataFileWriter1.write(al11);
                                    }

                                    written++;
                                    varincluded = true;
                                }

                            }
                        }
                    }

                    String incvar = "NotWritten";
                    if (varincluded) {

                        incvar = "Written";

                    }
                    String log = var.toString() + "\t" + Strings.concat(var.getAlleles(), Strings.semicolon) + "\t" + var.isIndel() + "\t" + var.isMultiallelic() + "\t" + var.getMAF() + "\t" + missingness + "\t" + missing + "\t" + sumdepth + "\t" + incvar;
                    logout.writeln(log);
                }
                if (ctr % 10000 == 0) {
                    //
                    double perc = (double) written / ctr;
                    perc *= 100;
                    DecimalFormat format = new DecimalFormat("##.##");
                    String strperc = format.format(perc);
                    System.out.println(
                            "File: " + f + "\t" +
                                    ctr + "\tlines, " +
                                    nrmultiallelic + " multi, " +
                                    failfilter + " filter, " +
                                    failinbreeding + " inbred, " +
                                    nrlowimpqual + " low R2, " +
                                    variantswithmissingness + " miss, " +
                                    lowmaf + " MAF, " +
                                    nrindels + " indels, " +
                                    written + " written (" + strperc + "%) - " +
                                    vcffile);
                }
                ctr++;
                ln = tf.readLine();
            }
            tf.close();
        }


        logout.close();

        System.out.println();
        System.out.println("Done");
        genotypeDataFileWriter1.close();


        snpFileWriter1.close();
        snpMapFileWriter.close();
        System.out.println(written + " written, " + nrvarswithdosage + " with dosages. " + lowmaf + " maf<" + mafthreshold + ", " + nrindels + " indels ," + nrmultiallelic + " multiallelic, " + variantswithmissingness + " with MAF>" + mafthreshold + " and high missigness " + vcffile);
    }
}
