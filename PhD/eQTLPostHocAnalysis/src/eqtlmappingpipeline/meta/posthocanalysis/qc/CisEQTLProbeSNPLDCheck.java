/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.qc;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.SortableSNP;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CisEQTLProbeSNPLDCheck {

    String probeSnpAnnotationfileLocation = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-SNPsPerProbe.txt";
    String eqtlfileLocation = "/Volumes/iSnackHD/Data/Projects/Isis/eQTLResults/2012-10-30-eQTLs-Cis/SNP-Initial/eQTLProbesFDR0.05-MetaIDs.txt";
    String qcoutdir = "/Volumes/iSnackHD/Data/Projects/Isis/eQTLResults/2012-10-30-eQTLs-Cis/SNP-Initial/QC2/";
    String kgLocation = "/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/";
    private String haplocation = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";

    
    double cr = 0.95;
    double hwe = 0.001;
    double maf = 0.01;
    public static void main(String[] args) {
        try {

            CisEQTLProbeSNPLDCheck p = new CisEQTLProbeSNPLDCheck();
            String outdir = "/Volumes/iSnackHD/tmp/";
            String hap = "/Data/1000g/";
            String snpp = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            p.determineSNPProbePairsWhichMayHaveFalsePositiveEffect(snpp, hap, outdir);
//            p.run();
//            p.calculateOverlap();
//	    p.countSNPsWithoutAnnotation();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void determineSNPProbePairsWhichMayHaveFalsePositiveEffect(String probeTranslation, String reference, String outdir) throws IOException {

        ProbeTranslation pb = new ProbeTranslation();
        pb.load(probeTranslation);


        TriTyperGenotypeData ds = new TriTyperGenotypeData();
        ds.load(reference);

        String[] snps = ds.getSNPs();

        TextFile output = new TextFile(outdir + "SNPProbeCombosWithPossibleHybArtifacts1Kg.txt", TextFile.W);
        for (byte chr = 1; chr < 23; chr++) {
            ArrayList<SortableSNP> snpsOnChr = new ArrayList<SortableSNP>();
            ArrayList<Integer> probesOnChr = new ArrayList<Integer>();

            int nrProbes = pb.getNumProbes();
            for (int p = 0; p < nrProbes; p++) {
                if (pb.getProbeChr(p) == chr) {
                    probesOnChr.add(p);
                }
            }

            for (int s = 0; s < snps.length; s++) {
                if (ds.getChr(s) == chr) {
                    snpsOnChr.add(new SortableSNP(snps[s], s, chr, ds.getChrPos(s), SortableSNP.SORTBY.ID));
                }
            }

            System.out.println(snpsOnChr.size() + "\tSNPs on Chr " + chr);
            System.out.println(probesOnChr.size() + "\tProbes on Chr " + chr);
            Collections.sort(snpsOnChr);

            HashMap<Integer, HashSet<Integer>> snpsInProbes = new HashMap<Integer, HashSet<Integer>>();

            for (int p = 0; p < probesOnChr.size(); p++) {
                String actualAnnotation = pb.getActualMappingPosition(probesOnChr.get(p));
                String[] elems = actualAnnotation.split(":");

                for (String pos : elems) {
                    String[] elems2 = pos.split("-");
                    Integer start = Integer.parseInt(elems2[0]);
                    Integer stop = Integer.parseInt(elems2[1]);

                    for (SortableSNP s : snpsOnChr) {
                        int chrPos = s.chrpos;
                        if (chrPos >= start && chrPos <= stop) {
                            HashSet<Integer> snpsInProbe = snpsInProbes.get(probesOnChr.get(p));
                            if (snpsInProbe == null) {
                                snpsInProbe = new HashSet<Integer>();
                            }
                            snpsInProbe.add(s.id);
                            snpsInProbes.put(probesOnChr.get(p), snpsInProbe);
                        }
                    }
                }
            }
            System.out.println(snpsInProbes.size() + "\tprobes with SNPs");

            SNPLoader loader = ds.createSNPLoader();
            DetermineLD ldcalc = new DetermineLD();

            ProgressBar progress = new ProgressBar(snpsOnChr.size(), "Testing chr: " + chr);

            HashSet<Integer> snpsNotPassingQC = new HashSet<Integer>();
            for (SortableSNP s : snpsOnChr) {
                // for each probe within 1Mb, check
                if (!snpsNotPassingQC.contains(s.id)) {
                    SNP snpObj1 = ds.getSNPObject(s.id);
                    loader.loadGenotypes(snpObj1);

                    if (snpObj1.getMAF() > maf && snpObj1.getHWEP() > hwe && snpObj1.getCR() > cr) {
                        for (int p = 0; p < probesOnChr.size(); p++) {
                            HashSet<Integer> snpsInProbe = snpsInProbes.get(probesOnChr.get(p));
                            if (snpsInProbe != null) {
                                boolean failsQC = false;
                                String qcStr = "";
                                for (Integer snpInProbe : snpsInProbe) {
                                    if (!snpsNotPassingQC.contains(snpInProbe)) {
                                        SNP snpObj2 = ds.getSNPObject(snpInProbe);
                                        loader.loadGenotypes(snpObj2);
                                        if (snpObj2.getMAF() > maf && snpObj2.getHWEP() > hwe && snpObj2.getCR() > cr) {
                                            Pair<Double, Double> ld = ldcalc.getLD(snpObj1, snpObj2, ds, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                                            if (ld.getLeft() > 0.2 || ld.getRight() > 0.2) {
                                                failsQC = true;
                                                qcStr += "\t" + snpObj2.getName() + " (" + ld.getLeft() + ", " + ld.getRight() + ")";
                                            }
                                        } else {
                                            snpsNotPassingQC.add(snpInProbe);
                                        }
                                        snpObj2.clearGenotypes();
                                    }
                                }
                                if (failsQC) {
                                    String outputStr = s.name + "\t" + pb.getProbes()[probesOnChr.get(p)] + qcStr;
                                    output.writeln(outputStr);
                                }
                            }
                        }
                    } else {
                        snpsNotPassingQC.add(s.id);
                    }
                    snpObj1.clearGenotypes();
                }
                progress.iterate();
            }
            progress.close();
        }
        output.close();
    }

    public void run() throws IOException {
        // load probes with snps
        Gpio.createDir(qcoutdir);

        TextFile probeSNPAnnotation = new TextFile(probeSnpAnnotationfileLocation, TextFile.R);
        String[] probeelems = probeSNPAnnotation.readLineElemsReturnObjects(TextFile.tab);

        HashMap<String, String[]> allProbesWithSNPs = new HashMap<String, String[]>();
        HashMap<String, Byte> probeToChrMap = new HashMap<String, Byte>();
        while (probeelems != null) {

            String probe = probeelems[0];
            Byte chr = ChrAnnotation.parseChr(probeelems[2]);
            if (probeelems.length > 8) {

                if (probeelems[8].trim().length() > 0) {
                    if (probeelems[8].startsWith(",")) {
                        probeelems[8] = probeelems[8].substring(1);
                    }
                    // probe has snps
                    String[] snps = Strings.comma.split(probeelems[8]);
                    allProbesWithSNPs.put(probe, snps);
                }
            }

            probeToChrMap.put(probe, chr);
            probeelems = probeSNPAnnotation.readLineElemsReturnObjects(TextFile.tab);
        }
        probeSNPAnnotation.close();

        System.out.println(allProbesWithSNPs.size() + " probes with SNPs.");

        DetermineLD ldcalc = new DetermineLD();
        HashSet<String> testedProbes = new HashSet<String>();
        int outctr = 0;

        HashMap<Byte, ArrayList<String[]>> eQTLsPerChromosome = new HashMap<Byte, ArrayList<String[]>>();
        TextFile eqtlfile = new TextFile(eqtlfileLocation, TextFile.R);
        String eqtlfileheader = eqtlfile.readLine();
        String[] eqtlfileelems = eqtlfile.readLineElems(TextFile.tab);
        while (eqtlfileelems != null) {

            Byte chr = ChrAnnotation.parseChr(eqtlfileelems[2]);
            ArrayList<String[]> e = eQTLsPerChromosome.get(chr);
            if (e == null) {
                e = new ArrayList<String[]>();
            }
            e.add(eqtlfileelems);


            eQTLsPerChromosome.put(chr, e);

            eqtlfileelems = eqtlfile.readLineElems(TextFile.tab);
        }
        eqtlfile.close();


        TextFile tpout = new TextFile(qcoutdir + "/TruePositives.txt", TextFile.W);
        tpout.writeln(eqtlfileheader);
        TextFile fpout = new TextFile(qcoutdir + "/FalsePositives.txt", TextFile.W);
        fpout.writeln(eqtlfileheader);
        TextFile upout = new TextFile(qcoutdir + "/UnknownPositives.txt", TextFile.W);

        TextFile log = new TextFile(qcoutdir + "/Log.txt", TextFile.W);
        log.writeln(eqtlfileheader + "\tSNPProbeDistance\tr2Sum\tmafsum\tSNPsIn1000G\tr2s");
        int nreffectsparsed = 0;
        for (byte chr = 1; chr < 23; chr++) {
//	for (byte chr = 2; chr < 3; chr++) {
            System.out.println("Chr " + chr);
            TriTyperGenotypeData kg = new TriTyperGenotypeData();
            kg.load(kgLocation + "/Chr" + chr + "/");

            TriTyperGenotypeData hap = new TriTyperGenotypeData();
            hap.load(haplocation);

            String[] snpsInHap = hap.getSNPs();
            ArrayList<SortableSNP> sortedSNPsInHapForChr = new ArrayList<SortableSNP>();
            for (int s = 0; s < snpsInHap.length; s++) {
                if (hap.getChr(s) == chr) {
                    sortedSNPsInHapForChr.add(new SortableSNP(snpsInHap[s], s, hap.getChr(s), hap.getChrPos(s), SortableSNP.SORTBY.CHRPOS));
                }
            }

            System.out.println(sortedSNPsInHapForChr.size() + " SNPs in HAP for Chr " + chr);
            Collections.sort(sortedSNPsInHapForChr);

            SNPLoader loader = kg.createSNPLoader();
            SNPLoader haploader = hap.createSNPLoader();

//	    TextFile eqtlfile = new TextFile(eqtlfileLocation, TextFile.R);
//	    elems = eqtlfile.readLineElemsReturnObjects(TextFile.tab); // header
//	    elems = eqtlfile.readLineElemsReturnObjects(TextFile.tab);
            int lnctr = 1;

            ArrayList<String[]> eqtls = eQTLsPerChromosome.get(chr);

            if (eqtls != null) {
                for (String[] elems : eqtls) {
                    String eSNP = elems[1];
                    Integer eSNPPos = Integer.parseInt(elems[3]);
                    String probe = elems[4];

                    if (probeToChrMap.get(probe) != null && !probeToChrMap.get(probe).equals(chr)) {
//		    System.out.println("Error: probe should map to chr: "+chr+" but it doesnt: "+probeToChrMap.get(probe)+"\t"+Strings.concat(elems, Strings.space));
                        upout.writeln(Strings.concat(elems, Strings.tab) + "\t" + "Probe has no mapping in new annotation");

                    } else if (probeToChrMap.get(probe) != null && probeToChrMap.get(probe).equals(chr)) {
//		    probeSNPsAlreadyOutput.add(eSNP+"-"+probe);
                        if (!allProbesWithSNPs.containsKey(probe)) {
                            tpout.writeln(Strings.concat(elems, Strings.tab));
                            log.writeln(elems[0] + "\t" + elems[1] + "\t" + elems[4] + "\tNo SNPs in Probe");
                        } else {

                            // probe contains SNPs
                            boolean probePassesQC = true;

                            Integer probePos = Integer.parseInt(elems[6]);

                            String[] probeSNPs = allProbesWithSNPs.get(probe);


                            Integer eSNPObjId = kg.getSnpToSNPId().get(eSNP);

                            String snpsInProbe = Strings.concat(probeSNPs, Strings.comma);

                            String qcString = "" + probeSNPs.length + "\t" + snpsInProbe;
                            String[] r2s = new String[probeSNPs.length];

                            double mafsum = 0;

                            double r2sum = 0;
                            int r2calc = 0;
                            if (eSNPObjId == null) {


                                // eSNP is not in 1000genomes.. take one up and one downstream, and see if they are in perfect LD, then take either one of them.

                                Integer hapeSNPId = hap.getSnpToSNPId().get(eSNP);
                                if (hapeSNPId == null) {
                                    log.writeln(eSNP + "\t not present in hapmap???");
                                    System.out.println(eSNP + "\t not present in hapmap???");
                                } else {
                                    SNP hapeSNP = hap.getSNPObject(hapeSNPId);
                                    haploader.loadGenotypes(hapeSNP);

                                    Integer selectedSNP = null;
                                    for (int i = 0; i < sortedSNPsInHapForChr.size(); i++) {
                                        if (eSNP.equals(sortedSNPsInHapForChr.get(i).name)) {
                                            selectedSNP = i;
                                            break;
                                        }
                                    }

                                    if (selectedSNP == null) {
                                        System.out.println(sortedSNPsInHapForChr.size());
                                        System.out.println(eSNP + " has no location in hapmap for chr " + chr);
                                    }

                                    boolean printoutput = false;

                                    // now go upstream and downstream until the r2 is no longer 1.

                                    Integer upstream = 0;
                                    Integer downstream = 0;

                                    // first look downstream...
                                    for (int i = selectedSNP + 1; (i < selectedSNP + 50 && i < sortedSNPsInHapForChr.size()); i++) {

                                        if (selectedSNP + 50 < sortedSNPsInHapForChr.size()) {
                                            String hapSNP2Name = sortedSNPsInHapForChr.get(i).name;

                                            if (kg.getSnpToSNPId().get(hapSNP2Name) != null) {
                                                SNP hapeSNP2 = hap.getSNPObject(sortedSNPsInHapForChr.get(i).id);
                                                haploader.loadGenotypes(hapeSNP2);

                                                // calculate LD. If it is 1, we have a candidate.
                                                double r2 = ldcalc.getRSquared(hapeSNP, hapeSNP2, hap, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                                                if (printoutput) {
                                                    System.out.println(selectedSNP + "\t" + sortedSNPsInHapForChr.get(selectedSNP).chrpos + "\t" + i + "\t" + sortedSNPsInHapForChr.get(i).chrpos + "\t" + eSNP + "\t" + hapSNP2Name + "\t" + r2 + "\t" + hapeSNP2.getMAF());
                                                }
                                                hapeSNP2.clearGenotypes();
                                                if (r2 >= 1d && eSNPObjId == null) {
                                                    eSNPObjId = kg.getSnpToSNPId().get(hapSNP2Name);
                                                    break;
                                                }
//					    else if (r2 <= 0.95) {
//						break;
//					    }
                                            }
                                        }
                                    }

                                    if (printoutput) {
                                        System.out.println("---");
                                    }
                                    // now check upstream, if we don't yet have a candidate...
                                    if (eSNPObjId == null) {
                                        for (int i = selectedSNP - 1; (i > selectedSNP - 50 && i > -1); i--) {
                                            String hapSNP2Name = sortedSNPsInHapForChr.get(i).name;

                                            if (kg.getSnpToSNPId().get(hapSNP2Name) != null) {
                                                SNP hapeSNP2 = hap.getSNPObject(sortedSNPsInHapForChr.get(i).id);
                                                haploader.loadGenotypes(hapeSNP2);

                                                // calculate LD. If it is 1, we have a candidate.
                                                double r2 = ldcalc.getRSquared(hapeSNP, hapeSNP2, hap, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                                                if (printoutput) {
                                                    System.out.println(selectedSNP + "\t" + sortedSNPsInHapForChr.get(selectedSNP).chrpos + "\t" + i + "\t" + sortedSNPsInHapForChr.get(i).chrpos + "\t" + eSNP + "\t" + hapSNP2Name + "\t" + r2 + "\t" + hapeSNP2.getMAF());
                                                }
                                                hapeSNP2.clearGenotypes();
                                                if (r2 >= 1d && eSNPObjId == null) {
                                                    eSNPObjId = kg.getSnpToSNPId().get(hapSNP2Name);
                                                    break;
                                                }
//					    else if (r2 <= 0.95) {
//						break;
//					    }
                                            }

                                        }
                                    }


//                                if (printoutput) {
//                                    System.out.println("");
//                                    System.out.println("");
//                                }

                                    if (eSNPObjId == null) {
                                        System.out.println("Could not find a proper proxy for\t" + eSNP);




                                        log.writeln(elems[0] + "\t" + elems[1] + "\t" + elems[4] + "\tCould not find a proper proxy for eSNP\t" + eSNP);
                                    } else {
                                        System.out.println("Found proxy for SNP " + eSNP + "\t" + kg.getSNPs()[eSNPObjId]);
                                    }
//                                if (printoutput) {
//                                    Integer id = hap.getSnpToSNPId().get("rs2240538");
//                                    SNP hapeSNP2 = hap.getSNPObject(id);
//                                    haploader.loadGenotypes(hapeSNP2);
//
//                                    double r2 = ldcalc.getRSquared(hapeSNP, hapeSNP2, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
//
//                                    Integer snppos = null;
//                                    for (int i = 0; i < sortedSNPsInHapForChr.size(); i++) {
//                                        if ("rs2240538".equals(sortedSNPsInHapForChr.get(i).name)) {
//                                            snppos = i;
//                                            break;
//                                        }
//                                    }
//
//                                    System.out.println(selectedSNP + "\t" + snppos + "\t" + eSNP + "\t" + hapeSNP2.getName() + "\t" + r2);
//
////				    System.exit(0);
//                                }
                                    hapeSNP.clearGenotypes();
//			    System.out.println("eSNP not present in 1kg " + eSNP);
                                }
                            }

                            if (eSNPObjId == null) {
                                upout.writeln(Strings.concat(elems, Strings.tab) + "\t" + "eSNP not present in 1kg and could not find a suitable proxy in HapMap2:" + eSNP);
                                log.writeln(elems[0] + "\t" + elems[1] + "\t" + elems[4] + "\tCould not find a proper proxy for eSNP\t" + eSNP + "\tcalculating LD in HapMap");
                                outctr++;
                            } else {
                                SNP eSNPObj = kg.getSNPObject(eSNPObjId);

                                loader.loadGenotypes(eSNPObj);

                                int snpctr = 0;
                                for (String probeSNP : probeSNPs) {
                                    Integer pSNPObjId = kg.getSnpToSNPId().get(probeSNP);
                                    if (pSNPObjId == null) {
//				    System.out.println("pSNP not present in 1kg " + probeSNP);
                                        r2s[snpctr] = "-";
                                    } else {
                                        SNP pSNPObj = kg.getSNPObject(pSNPObjId);

                                        if (pSNPObj == eSNPObj) {
//					System.out.println("FALSE POSITIVE:\tSNPs in perfect LD for probe " + probe);
                                            r2s[snpctr] = "1.0";
                                            probePassesQC = false;
                                        } else {

                                            loader.loadGenotypes(pSNPObj);

                                            double r2 = ldcalc.getRSquared(eSNPObj, pSNPObj, hap, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);

                                            if (r2 > 0.2) {
                                                probePassesQC = false;
//					    System.out.println("FALSE POSITIVE:\t" + eSNPObj.getGenotypes().length + "\t" + eSNP + "\t" + pSNP + "\t" + probe + "\t" + r2);
                                            }
                                            mafsum += pSNPObj.getMAF();

                                            r2sum += r2;
                                            r2calc++;
                                            r2s[snpctr] = "" + r2;

                                            pSNPObj.clearGenotypes();

                                        }
                                    }
                                    snpctr++;
                                }

                                int distance = eSNPObj.getChrPos() - probePos;

                                if (r2calc > 0) {
                                    r2sum /= r2calc;
                                    mafsum /= r2calc;
                                }
//String s = "SNPProbeDistance\tr2Sum\tmafsum\tSNPsIn1000G\tr2s";
                                qcString = distance + "\t" + r2sum + "\t" + mafsum + "\t" + qcString + "\t" + Strings.concat(r2s, Strings.comma);

                                if (!probePassesQC) {
                                    fpout.writeln(Strings.concat(elems, Strings.tab));
                                    log.writeln(elems[0] + "\t" + elems[1] + "\t" + elems[4] + "\t" + qcString);
                                    outctr++;
                                } else {
                                    tpout.writeln(Strings.concat(elems, Strings.tab));
                                    log.writeln(elems[0] + "\t" + elems[1] + "\t" + elems[4] + "\t" + qcString);
                                    outctr++;
                                }
                                eSNPObj.clearGenotypes();
                            }
                            testedProbes.add(probe);
                        }
                    } else {
                        outctr++;
                        upout.writeln(Strings.concat(elems, Strings.tab) + "\t" + "Probe has null chromosome");
                        log.writeln(elems[0] + "\t" + elems[1] + "\t" + elems[4] + "\t" + "Probe has null chromosome");
//		    System.out.println("Chr for probe "+probe+" is null: "+probeToChrMap.get(probe));
//		    System.out.println(Strings.concat(elems, Strings.space));
//		} else {
////		    upout.writeln(Strings.concat(elems, Strings.tab) + "\t" + "Probe has null chromosome");
                    }
                }
            }


            haploader.close();
            loader.close();
        }

        System.out.println(outctr + "\teQTLs tested");

        log.close();
        tpout.close();
        fpout.close();
        upout.close();
    }

    private void calculateOverlap() throws IOException {
        // determine for each file the number of uniquely affected probes and snps/...

        HashSet<String> uniqueSNPsInTruePositives = new HashSet<String>();
        HashSet<String> uniqueSNPsInFalsePositives = new HashSet<String>();
        HashSet<String> uniqueSNPsInUnknownPositives = new HashSet<String>();


        HashSet<String> uniqueProbesInTruePositives = new HashSet<String>();
        HashSet<String> uniqueProbesInFalsePositives = new HashSet<String>();
        HashSet<String> uniqueProbesInUnknownPositives = new HashSet<String>();

        TextFile tf = new TextFile(qcoutdir + "TruePositives.txt", TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            String probe = elems[4];
            String snp = elems[1];

            uniqueSNPsInTruePositives.add(probe);
            uniqueProbesInTruePositives.add(snp);

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        tf = new TextFile(qcoutdir + "FalsePositives.txt", TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            String probe = elems[4];
            String snp = elems[1];

            uniqueSNPsInFalsePositives.add(probe);
            uniqueProbesInFalsePositives.add(snp);


            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        tf = new TextFile(qcoutdir + "UnknownPositives.txt", TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            String probe = elems[4];
            String snp = elems[1];

            uniqueSNPsInUnknownPositives.add(probe);
            uniqueProbesInUnknownPositives.add(snp);

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        System.out.println(uniqueSNPsInTruePositives.size() + "\tunique probes in true positives");
        System.out.println(uniqueSNPsInFalsePositives.size() + "\tunique probes in false positives");
        System.out.println(uniqueSNPsInUnknownPositives.size() + "\tunique probes in unknown positives");


        System.out.println(uniqueProbesInTruePositives.size() + "\tunique probes in true positives");
        System.out.println(uniqueProbesInFalsePositives.size() + "\tunique probes in false positives");
        System.out.println(uniqueProbesInUnknownPositives.size() + "\tunique probes in unknown positives");



        // now determine the overlap.....

        String[] probes = uniqueProbesInTruePositives.toArray(new String[0]);
        int ctr = 0;
        for (int i = 0; i < probes.length; i++) {
            if (uniqueProbesInFalsePositives.contains(probes[i])) {
                ctr++;
            }
        }

        System.out.println(ctr + "\tprobes overlap between true positives and false positives");
        String[] snps = uniqueSNPsInTruePositives.toArray(new String[0]);
        ctr = 0;
        for (int i = 0; i < snps.length; i++) {
            if (uniqueSNPsInFalsePositives.contains(snps[i])) {
                ctr++;
            }
        }
        System.out.println(ctr + "\tsnps overlap between true positives and false positives");
    }

    public void countSNPsWithoutAnnotation() throws IOException {
        int ctr = 0;
        int total = 0;

        for (byte chr = 1; chr < 23; chr++) {
//	for (byte chr = 2; chr < 3; chr++) {
            System.out.println("Chr " + chr);
            TriTyperGenotypeData kg = new TriTyperGenotypeData();
            kg.load(kgLocation + "/Chr" + chr + "/");
            String[] snps = kg.getSNPs();

            for (int i = 0; i < snps.length; i++) {
                if (kg.getChr(i) != chr) {
//		    System.out.println(snps[i] + "\tshould belong to chr " + chr + " but has no such annotation");
                    ctr++;
                }


            }
            total += snps.length;
        }

        TextFile out = new TextFile("/Volumes/BackupDisk/Liftover/snps.bed", TextFile.W);
        for (byte chr = 1; chr < 23; chr++) {
            TextFile in = new TextFile(kgLocation + "/Chr" + chr + "/SNPMappings_b37.txt", TextFile.R);
            String[] elems = in.readLineElems(TextFile.tab);
            while (elems != null) {
                out.writeln("Chr" + elems[0] + " " + elems[1] + " " + elems[1] + " " + elems[2]);
            }
            in.close();
        }

        out.close();
        System.out.println(ctr + " out of " + total + " have no annotation...");
    }
}
