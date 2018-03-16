/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.probes;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Chromosome;
import umcg.genetica.containers.Exon;
import umcg.genetica.containers.Gene;
import umcg.genetica.containers.Transcript;
import umcg.genetica.ensembl.Features;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ProbeFilterAfterMapping {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        // load sequence matched to probe

        System.out.println("Running probe filter..");
        try {
            ProbeFilterAfterMapping p = new ProbeFilterAfterMapping();
            System.out.println("Reading old annotation file for sequences");
            // 
//            p.readSequencesFromOldAnnotationFile("D:\\Skydrive\\MetaAnalysisAnnotationFiles\\2011-10-06-ProbeTranslationTable+H8HT12Conversion.log_reannotatedHG18_96PercIdentity.txt");
            p.readSequencesFromOldAnnotationFile("/Volumes/iSnackHD/Skydrive/MetaAnalysisAnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotationFormat.txt");
//
//	    // convert probe mapping file to new probe translation format
//	    p.convertProbeMappingFileToProbeAnnotationFile("/Data/ProbeAnnotation/ProbeMappings/2012-04-23-IlluminaAll96PercentIdentity/all_uniqelymapping",
//		    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotationFormat.txt");
//
//	    // filter out probes which have < 50 mapping length
//	    System.out.println("Removing probes with < 48 mapping");
//	    p.removeProbesWithIncorrectMappingLength("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotationFormat.txt",
//		    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut.txt", 48, 50);


            // determine where the non-mapping probes come from: are they mostly exon-junction spanning, or hypothetical probes
//	    p.annotateNonMappingProbesFromIlluminaFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-02-22-ProbeAnnotation96PercIdentity-ExonJunctionMappingsFixed-ProbesWithWrongMappingLengthFilteredOut.txt",
//		    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-02-22-ProbeAnnotation96PercIdentity-ExonJunctionMappingsFixed.txt",
//		    "/Data/ProbeAnnotation/ProbeAnnotationAccordingToIllumina/Illumina/");

            // annotate probes for genes, snps, TSS, etc
//2012-03-07-ProbeAnnotation96PercIdentity-ExonJunctionMappingsFixed-TranscriptMappingsFixed-ProbesWithWrongMappingLengthFilteredOut.txt
//            p.loadSNPAnnotation("/Data/SNPReferenceData/1000G-20110521-TriTyper/Chr1/SNPMappings.txt");
//            p.annotateProbesWithSNPData("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut.txt",
//                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-SNPsPerProbe.txt");
//	    // output final probe annotation
//
//
            p.annotateProbesWithEnsemblGenes("/Volumes/iSnackHD/Skydrive/MetaAnalysisAnnotationFiles/structures_b54.txt",
                    "/Volumes/iSnackHD/Skydrive/MetaAnalysisAnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut.txt",
                    "/Volumes/iSnackHD/Skydrive/MetaAnalysisAnnotationFiles/2012-08-08-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation.txt");
//
//	    p.convertToProbeTranslationFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-14-ProbeTranslationTable+H8HT12Conversion.txt",
//		    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation.txt",
//		    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt");

//	    p.convertToPlatformAnnotationFiles("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-03-07-ProbeAnnotation96PercIdentity-ExonJunctionMappingsFixed-TranscriptMappingsFixed-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt",
//		    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    private HashMap<String, String> oldSeq;
    private ArrayList<String> probeNrs;
    HashMap<Byte, HashSet<Integer>> chromosomePositionsWithSNPs = new HashMap<Byte, HashSet<Integer>>();
    HashMap<Byte, HashMap<Integer, String>> chromosomeSNPPositions = new HashMap<Byte, HashMap<Integer, String>>();

    private void readSequencesFromOldAnnotationFile(String volumesData2MarjoleinHomeAccountmarjol) throws IOException {
        oldSeq = new HashMap<String, String>();

        probeNrs = new ArrayList<String>();
        TextFile old = new TextFile(volumesData2MarjoleinHomeAccountmarjol, TextFile.R);

        String[] elems = old.readLineElems(TextFile.tab);
        while (elems != null) {
            oldSeq.put(elems[0], elems[1]);
            probeNrs.add(elems[0]);
            elems = old.readLineElems(TextFile.tab);
        }

        old.close();
    }

    public void convertProbeMappingFileToProbeAnnotationFile(String in, String out) throws IOException {

        TextFile tf = new TextFile(in, TextFile.R);

        TextFile tfout = new TextFile(out, TextFile.W);


        String[] elems = tf.readLineElems(TextFile.tab);

        HashMap<String, String> annotation = new HashMap<String, String>();

        while (elems != null) {

            if (elems.length > 5) {
                String probe = elems[1].trim();
                String chr = elems[3].trim();
                String chrstart = elems[4].trim();
                String chrend = elems[5].trim();

                String[] chrStartElems = chrstart.split(",");
                String[] chrEndElems = chrend.split(",");


                String finalAnnot = "";

                for (int i = 0; i < chrEndElems.length; i++) {
                    if (i == 0) {
                        finalAnnot += chrStartElems[i].trim() + "-" + chrEndElems[i].trim();
                    } else {
                        finalAnnot += ":" + chrStartElems[i].trim() + "-" + chrEndElems[i].trim();
                    }
                }


                String outStr = chr + "\t" + finalAnnot + "\t-";

                System.out.println("Parsed: " + probe + "\t" + chr + "\t" + finalAnnot);

                annotation.put(probe, outStr);
            } else {
                String ln = Strings.concat(elems, Strings.tab);
                System.out.println("Could not parse: " + ln);
            }


            elems = tf.readLineElems(TextFile.tab);
        }


        for (int i = 0; i < probeNrs.size(); i++) {
            String probe = probeNrs.get(i).trim();
            String seq = oldSeq.get(probe).trim();
            String annot = annotation.get(probe);
            if (annot == null) {
                annot = "-\t-\t-";
//		System.out.println("No annotation for probe: "+probe);
            } else {
            }
            String outStr = probe + "\t" + seq + "\t" + annot;
            System.out.println(outStr);
            tfout.writeln(outStr);
        }



        tfout.close();
        tf.close();

    }

    private void removeProbesWithIncorrectMappingLength(String probetranslationfile, String out, int minlength, int maxlength) throws IOException {


        ProbeTranslation p = new ProbeTranslation();
        p.load(probetranslationfile);
        TextFile outfile = new TextFile(out, TextFile.W);

        String[] probes = p.getProbes();

        System.out.println("Probe\tseq\tchr\tchrpos\thugo");
        outfile.writeln("Probe\tseq\tchr\tchrpos\thugo");
        int probesWithMapping = 0;

        for (int i = 0; i < probes.length; i++) {

            byte chr = p.getProbeChr(i);
            int snpcounter = 0;
            int probelengthsum = 0;
            ArrayList<String> snpnamesMappingToProbe = new ArrayList<String>();
            int midpos = 0;
            int snpvicinitycounter = 0;
            String chrpos = p.getActualMappingPosition(i);
            if (chr > 0) {
                String[] list = chrpos.split(":"); // start1-end1:start2-end2
                try {

                    for (String s : list) {
                        String[] list2 = s.split("-");
                        if (list2.length == 2) {
                            Integer start = Integer.parseInt(list2[0]);
                            Integer end = Integer.parseInt(list2[1]);
                            int diff = end - start + 1;
                            probelengthsum += diff;
                            midpos += start + ((int) Math.abs((double) (start - end) / 2));
                        }

                    }

                    midpos /= list.length;

                } catch (Exception e2) {
                    // no mapping available
                }
            }

            String seq = oldSeq.get(probes[i]);
            if (probelengthsum >= minlength && probelengthsum <= maxlength) {
                probesWithMapping++;
//		System.out.println(probelengthsum + "\t" + probes[i] + "\t" + seq + "\t" + ChrAnnotation.parseByte(chr) + "\t" + chrpos);
                outfile.writeln(probes[i] + "\t" + seq + "\t" + ChrAnnotation.parseByte(chr) + "\t" + chrpos + "\t-");
            } else {
                System.out.println(probelengthsum + "\t" + probes[i] + "\t" + seq + "\t-\t-");
                outfile.writeln(probes[i] + "\t" + seq + "\t-\t-\t-");
            }
        }
        System.out.println(probesWithMapping + " out of " + probes.length + " probes have a mapping " + ((double) probesWithMapping / probes.length));
        outfile.close();
    }

    private void annotateNonMappingProbesFromIlluminaFile(String probeannotationfile, String myOriginalMappings, String illuminaannotationDir) throws IOException {

        ArrayList<String> sequencesOfProbesThatDontMap = new ArrayList<String>();
        HashMap<String, String> myMappings = new HashMap<String, String>();

        TextFile tf = new TextFile(probeannotationfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
//	    if (elems[2].equals("-")) {
            sequencesOfProbesThatDontMap.add(elems[1]);
//	    }

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        TextFile mym = new TextFile(myOriginalMappings, TextFile.R);

        elems = mym.readLineElems(TextFile.tab);

        while (elems != null) {
            myMappings.put(elems[1], (elems[2] + "-" + elems[3]));
            elems = mym.readLineElems(TextFile.tab);
        }

        mym.close();

        HashSet<String> seqs = new HashSet<String>();
        seqs.addAll(sequencesOfProbesThatDontMap);

        String[] annotationfiles = Gpio.getListOfFiles(illuminaannotationDir);


        for (String f : annotationfiles) {
            TextFile annotationfile = new TextFile(illuminaannotationDir + f, TextFile.R);

            elems = annotationfile.readLineElems(TextFile.tab);

            int predicted = 0;
            int total = 0;
            int exonjunctionspanning = 0;
            int probesWithoutIlluminaMapping = 0;

            while (elems != null) {
                if (elems.length >= 17) {
                    if (seqs.contains(elems[17])) {
                        String mapping = myMappings.get(elems[17]);
                        if (elems[22].toLowerCase().contains("hypothetical") || elems[22].toLowerCase().contains("predicted")) {

                            predicted++;

                        } else if (elems[20].contains(":")) {
                            exonjunctionspanning++;
                        }

                        if (elems[20].trim().length() == 0) {
                            probesWithoutIlluminaMapping++;
                        }

//			System.out.println(mapping + "\t" + f + "\t" + elems[17] + "\t" + elems[18] + "\t" + elems[20] + "\t" + elems[4] + "\t" + elems[22]);


                        total++;
                    }
                }
                elems = annotationfile.readLineElems(TextFile.tab);
            }

            System.out.println(predicted + "\t" + total + "\t" + ((double) predicted / total));
            System.out.println(exonjunctionspanning);
            System.out.println(probesWithoutIlluminaMapping);
            System.out.println("");
            annotationfile.close();
        }


    }

    private void loadSNPAnnotation(String dbsnp) throws IOException {
        TextFile snps = new TextFile(dbsnp, TextFile.R);

        String[] elems = snps.readLineElems(TextFile.tab);

        while (elems != null) {

            String chr = elems[0];
            Integer pos = Integer.parseInt(elems[1]);
            byte snpChr = ChrAnnotation.parseChr(chr);
            HashSet<Integer> snpsOnChr = chromosomePositionsWithSNPs.get(snpChr);
            HashMap<Integer, String> snpNames = chromosomeSNPPositions.get(snpChr);
            if (snpsOnChr == null) {
                snpNames = new HashMap<Integer, String>();
                snpsOnChr = new HashSet<Integer>();
            }
            snpsOnChr.add(pos);

            snpNames.put(pos, elems[2]);
            chromosomePositionsWithSNPs.put(snpChr, snpsOnChr);
            chromosomeSNPPositions.put(snpChr, snpNames);
            elems = snps.readLineElems(TextFile.tab);
        }

        snps.close();


    }

    private void annotateProbesWithSNPData(String input, String out) throws IOException {
        ProbeTranslation p = new ProbeTranslation();
        p.load(input);
        TextFile outfile = new TextFile(out, TextFile.W);
        String[] probes = p.getProbes();

        int probeswithsnps = 0;
        int probeswithannotation = 0;

        System.out.println("probe\tchr\tmidpos\tlengthofmapping\tnumsnps\tsnpswithin250kbofmidpos\tsnpsmapping");
        outfile.writeln("probe\tchr\tlengthofmapping\tnumsnps\tsnpswithin250kbofmidpos\tsnpsmapping");
        ProgressBar pb = new ProgressBar(probes.length);
        for (int i = 0; i < probes.length; i++) {

            byte chr = p.getProbeChr(i);
            int snpcounter = 0;
            int probelengthsum = 0;
            ArrayList<String> snpnamesMappingToProbe = new ArrayList<String>();
            int midpos = 0;
            int snpvicinitycounter = 0;
            String chrpos = p.getActualMappingPosition(i);
            if (chr > 0) {


//		System.out.println(i+"\t"+chr+"\t"+chrpos);

                String[] list = chrpos.split(":"); // start1-end1:start2-end2

                HashSet<Integer> snpsForChr = chromosomePositionsWithSNPs.get(chr);

                HashMap<Integer, String> snpsForChrNames = chromosomeSNPPositions.get(chr);



                try {


                    for (String s : list) {
                        String[] list2 = s.split("-");
                        if (list2.length == 2) {
                            Integer start = Integer.parseInt(list2[0]);
                            Integer end = Integer.parseInt(list2[1]);
                            int diff = end - start + 1;
                            probelengthsum += diff;
                            midpos += start + ((int) Math.abs((double) (start - end) / 2));
                            for (int pos = start; pos < end + 1; pos++) {
                                if (snpsForChr.contains(pos)) {
                                    snpcounter++;

                                    String name = snpsForChrNames.get(pos);
                                    snpnamesMappingToProbe.add(name);
                                }
                            }
                        }

                    }

                    midpos /= list.length;

                    for (int pos = midpos; pos < (midpos + 250000); pos++) {
                        if (snpsForChr.contains(pos)) {
                            snpvicinitycounter++;
                        }
                    }

                    for (int pos = midpos; pos > (midpos - 250000); pos--) {
                        if (snpsForChr.contains(pos)) {
                            snpvicinitycounter++;
                        }
                    }



                } catch (Exception e2) {
                    // no mapping available
                }


            }


            String allSNPNames = Strings.concat(snpnamesMappingToProbe, Strings.comma);

            if (snpcounter > 0) {
                probeswithsnps++;
            }

            outfile.writeln(probes[i] + "\t" + oldSeq.get(probes[i]) + "\t" + ChrAnnotation.parseByte(chr) + "\t" + chrpos + "\t" + midpos + "\t" + probelengthsum + "\t" + snpcounter + "\t" + snpvicinitycounter + "\t" + allSNPNames);

            if (probelengthsum == 50) {
                probeswithannotation++;
//		System.out.println(i + "\t" + chr + "\t" + chrpos + "\t" + midpos + "\t" + probelengthsum + "\t" + snpcounter + "\t" + snpvicinitycounter + "\t" + allSNPNames);
            } else {
//		System.out.println(i + "\t-\t-\t"+ midpos + "\t" + probelengthsum + "\t-\t-\t-");
            }


            pb.set(i);
        }

        pb.close();

        System.out.println(probeswithsnps + "\t out of " + probeswithannotation + " tested probes contain SNPs");
        outfile.close();
    }

    private void annotateProbesWithEnsemblGenes(String featurefile, String annotationfile, String outfile) throws IOException {
        Features f = new Features();
        f.loadAnnotation(featurefile);

        ProbeTranslation p = new ProbeTranslation();
        p.load(annotationfile);

        ChrAnnotation chran = new ChrAnnotation();

        TextFile out = new TextFile(outfile, TextFile.W);
        out.writeln("Probe\tSeq\tChr\tPos\tEnsemblGene\tHugo\tEnsemblGeneLoc\tEnsemblGeneStrand\tEnsemblTranscript\tEnsemblTranscriptLoc\tEnsemblExon\tEnsemblExonLoc");

        String defaultString = "\t-\t-\t-\t-\t-\t-\t-";
        HashMap<String, Integer> probetranslationtable = p.getProbeTranslationTable();
        boolean[] presentOnHT12v3 = new boolean[p.getNumProbes()];
        Set<String> oldProbeIds = probetranslationtable.keySet();
        int probectr = 0;
        for (String s : oldProbeIds) {
            if (s.contains("HumanHT-12_V3_0_R2_11283641_A.txt")) {
                Integer probeId = probetranslationtable.get(s);
                presentOnHT12v3[probeId] = true;
                probectr++;
            }
        }

        System.out.println(probectr + " probes on HT12v3");
        int numProbes = p.getNumProbes();

        String[] probes = p.getProbes();
        for (int i = 0; i < numProbes; i++) {
            byte probeChr = p.getProbeChr(i);
            int ctr = 0;
            if (probeChr <= 0 || probeChr >= 25) {
                // do some stuff
//                out.writeln(probes[i] + "\t" + oldSeq.get(probes[i]));
            } else if (probeChr > 0 && probeChr < 25) {
                int probeMidpointPosition = p.getProbeChrPos(i);
                if (probeMidpointPosition <= 0) {
                    // do default stuff
//                    out.writeln(probes[i] + "\t" + oldSeq.get(probes[i]));
                } else if (probeMidpointPosition > 0) {
                    Chromosome c = f.getChromosomeHash().get(ChrAnnotation.parseByte(probeChr));
                    if (c == null) {
                        System.out.println("Chromsome " + probeChr + " not found for probe " + i);
                    }
                    HashMap<String, Gene> genes = c.getGenesHash();

                    Set<Map.Entry<String, Gene>> availableGenes = genes.entrySet();
                    String ensemblGenes = "";
                    String strands = "";
                    String hgnc = "";

                    String TStart = "";
                    String TEnd = "";

                    HashSet<Gene> genesMappingToRegion = new HashSet<Gene>();
                    for (Entry<String, Gene> geneset : availableGenes) {

                        Gene gene = geneset.getValue();
                        int start = gene.getStart();
                        int end = gene.getEnd();

                        if (probeMidpointPosition >= start && probeMidpointPosition <= end) {
                            genesMappingToRegion.add(gene);
                            ctr++;
                        }
                    }

                    Gene[] genesArr = genesMappingToRegion.toArray(new Gene[0]);
                    String originalChrPos = p.getActualMappingPosition(i);

                    // if this has more than 1 element, shit has hit the fan. (exonboundary spanning probe)
                    String[] probeChrExonPositions = originalChrPos.split(":"); // format: exon1-exon1end:exon2-exon2end
                    boolean probeSpansExons = false;
                    if (probeChrExonPositions.length > 1) {
                        probeSpansExons = true;
                    }


                    HashSet<Transcript> matchedTranscripts = new HashSet<Transcript>();
                    HashMap<Transcript, ArrayList<Exon>> matchedExonsPerTranscript = new HashMap<Transcript, ArrayList<Exon>>();


                    // find out whether the probe maps to an exon within a gene transcript..
                    for (Gene g : genesArr) { // loop genes
                        HashMap<String, Transcript> transcripts = g.getTranscripts();
                        Set<Entry<String, Transcript>> transcriptset = transcripts.entrySet();


                        HashMap<Transcript, Integer> exonPerTranscriptCounter = new HashMap<Transcript, Integer>();
                        HashSet<Transcript> selectedTranscripts = new HashSet<Transcript>();

                        // first check whether the transcript is within this gene...
                        int genestart = g.getStart();
                        int geneend = g.getEnd();
                        if (probeMidpointPosition >= genestart && probeMidpointPosition <= geneend) {
                            // probe maps within gene region. does it map to any transcript of the gene?

                            for (Entry<String, Transcript> transcriptEntry : transcriptset) { // look into each transcript
                                Transcript t = transcriptEntry.getValue();
                                Exon[] exonsForTranscript = t.getExonsRanked();

                                int transcriptStart = t.getStart();
                                int transcriptEnd = t.getEnd();

                                if (probeMidpointPosition >= transcriptStart && probeMidpointPosition <= transcriptEnd) {
                                    // probe maps to this exon, but now check whether it is intronic or exonic

                                    for (String exonbound : probeChrExonPositions) { // look per exon
                                        String[] chrposelems = exonbound.split("-"); // this gives the actual start and stop positions

                                        int exonboundstart = Integer.parseInt(chrposelems[0]);
                                        int exonboundstop = Integer.parseInt(chrposelems[1]);

                                        for (Exon e : exonsForTranscript) {
                                            int exonStart = e.getStart();
                                            int exonStop = e.getEnd();
                                            if (exonboundstart >= exonStart && exonboundstop <= exonStop) {
                                                // mathedExons.add(e);
                                                Integer exoncounter = exonPerTranscriptCounter.get(t);
                                                if (exoncounter == null) {
                                                    exoncounter = 0;
                                                }
                                                exoncounter++;
                                                selectedTranscripts.add(t);
                                                exonPerTranscriptCounter.put(t, exoncounter);
                                                ArrayList<Exon> matchedExonsForTs = matchedExonsPerTranscript.get(t);
                                                if (matchedExonsForTs == null) {
                                                    matchedExonsForTs = new ArrayList<Exon>();
                                                }
                                                matchedExonsForTs.add(e);
                                                matchedExonsPerTranscript.put(t, matchedExonsForTs);
                                            }
                                        }
                                    }
                                }


                            }

                            // we now know which transcript maps to which probe
                            // the number of matched exons for each transcript should be equal to the number of exons in the annotation..
                            for (Transcript ts : selectedTranscripts) {
                                Integer exonsMatchedForTs = exonPerTranscriptCounter.get(ts);
                                if (exonsMatchedForTs.equals(probeChrExonPositions.length)) {
                                    matchedTranscripts.add(ts); // I guess this is a hit...
                                }
                            }
                        }


                    }

                    String finalOutput = "";
                    if (matchedTranscripts.isEmpty()) {
                        // apparently, the probe does not map to an exon, although it could map to an intron of a transcript.. (what does it then probe for??)
                        matchedTranscripts = new HashSet<Transcript>();
                        HashSet<Gene> matchedGenesWithoutTranscript = new HashSet<Gene>();
                        for (Gene g : genesArr) { // loop genes
                            HashMap<String, Transcript> transcripts = g.getTranscripts();
                            Set<Entry<String, Transcript>> transcriptset = transcripts.entrySet();

                            // first check whether the transcript is within this gene...
                            int genestart = g.getStart();
                            int geneend = g.getEnd();
                            if (probeMidpointPosition >= genestart && probeMidpointPosition <= geneend) {
                                // it maps to this gene
                                for (Entry<String, Transcript> transcriptEntry : transcriptset) { // look into each transcript
                                    Transcript t = transcriptEntry.getValue();

                                    int transcriptStart = t.getStart();
                                    int transcriptEnd = t.getEnd();

                                    if (probeMidpointPosition >= transcriptStart && probeMidpointPosition <= transcriptEnd) {
                                        // it maps to this transcript
                                        matchedTranscripts.add(t);
                                    }
                                }

                                if (matchedTranscripts.isEmpty()) {
                                    // probe maps to gene region but not to transcript.
                                    matchedGenesWithoutTranscript.add(g);
                                }
                            }
                        }

                        if (matchedTranscripts.isEmpty()) {
                            // probe does not map to any known transcript..
                            if (matchedGenesWithoutTranscript.isEmpty()) {
                                // also no matching gene
                            } else {
                                String name = "";
                                String annot = "";
                                String loc = "";
                                for (Gene g : matchedGenesWithoutTranscript) {
                                    if (name.length() > 0) {
                                        name += "," + g.getName();
                                        annot += "," + g.getAnnotation();
                                        loc += "," + g.getStart() + "-" + g.getEnd();
                                    } else {
                                        name += g.getName();
                                        annot += g.getAnnotation();
                                        loc += g.getStart() + "-" + g.getEnd();
                                    }
                                }

                                finalOutput += name + "\t" + annot + "\t" + loc;
                            }
                        } else {
                            // probe apparently maps into an intron..
                            HashSet<Gene> matchedGenes = new HashSet<Gene>();
                            for (Transcript t : matchedTranscripts) {
                                matchedGenes.add(t.getParentGene());
                            }

                            String name = "";
                            String annot = "";
                            String loc = "";
                            for (Gene g : matchedGenes) {
                                if (name.length() > 0) {
                                    name += "," + g.getName();
                                    annot += "," + g.getAnnotation();
                                    loc += "," + g.getStart() + "-" + g.getEnd();
                                } else {
                                    name += g.getName();
                                    annot += g.getAnnotation();
                                    loc += g.getStart() + "-" + g.getEnd();
                                }
                            }

                            finalOutput += name + "\t" + annot + "\t" + loc;

                            name = "";
                            loc = "";
                            for (Transcript t : matchedTranscripts) {
                                if (name.length() > 0) {
                                    name += "," + t.getName();
                                    loc += "," + t.getStart() + "-" + t.getEnd();
                                } else {
                                    name += t.getName();
                                    loc += t.getStart() + "-" + t.getEnd();
                                }
                            }

                            finalOutput += "\t" + name + "\t" + loc;

                        }

                    } else {
                        // yippee, possibly found a match!
                        // probe maps within an exon of a gene transcript
                        HashSet<Gene> matchedGenes = new HashSet<Gene>();
                        for (Transcript t : matchedTranscripts) {
                            matchedGenes.add(t.getParentGene());
                        }

                        String name = "";
                        String annot = "";
                        String strand = "";
                        String loc = "";
                        for (Gene g : matchedGenes) {
                            if (name.length() > 0) {
                                name += "," + g.getName();
                                annot += "," + g.getAnnotation();
                                strand += "," + g.getStrand();
                                loc += "," + g.getStart() + "-" + g.getEnd();
                            } else {
                                name += g.getName();
                                annot += g.getAnnotation();
                                loc += g.getStart() + "-" + g.getEnd();
                                strand += g.getStrand();
                            }
                        }

                        finalOutput += name + "\t" + annot + "\t" + loc + "\t" + strand;

                        name = "";
                        loc = "";
                        for (Transcript t : matchedTranscripts) {
                            if (name.length() > 0) {
                                name += "," + t.getName();
                                loc += "," + t.getStart() + "-" + t.getEnd();
                            } else {
                                name += t.getName();
                                loc += t.getStart() + "-" + t.getEnd();
                            }
                        }

                        finalOutput += "\t" + name + "\t" + loc;
                        name = "";
                        loc = "";
                        for (Transcript t : matchedTranscripts) {
                            ArrayList<Exon> matchedExons = matchedExonsPerTranscript.get(t);
                            String eName = "";
                            String eLoc = "";
                            for (Exon e : matchedExons) {
                                if (eName.length() > 0) {
                                    eName += "," + e.getName();
                                    eLoc += "," + e.getStart() + "-" + e.getEnd();
                                } else {
                                    eName += e.getName();
                                    eLoc += e.getStart() + "-" + e.getEnd();
                                }
                            }

                            if (name.length() > 0) {
                                name += ";" + eName;
                                loc += ";" + eLoc;
                            } else {
                                name += eName;
                                loc += eLoc;
                            }
                        }

                        finalOutput += "\t" + name + "\t" + loc;
                    }








                    if (finalOutput.length() > 0) {
                        out.writeln(probes[i] + "\t" + oldSeq.get(probes[i]) + "\t" + probeChr + "\t" + p.getActualMappingPosition(i) + "\t" + finalOutput);
                    } else {
                        out.writeln(probes[i] + "\t" + oldSeq.get(probes[i]) + "\t" + probeChr + "\t" + p.getActualMappingPosition(i) + defaultString);
                    }


                }
            }
        }
        out.close();
    }

// we have now found the set of transcripts to which the probe maps.. now determine which  
//                                if (mathedExons.size() == probeChrExonPositions.length) {
//                                    // we found a match :) now determine the transcripts that are shared by this set of exons
//                                    Exon[] matchedExonArr = mathedExons.toArray(new Exon[0]);
//
//                                    HashMap<Transcript, Integer> transcriptCounter = new HashMap<Transcript, Integer>();
//                                    HashSet<Transcript> allTranscripts = new HashSet<Transcript>();
//                                    for (Exon e : matchedExonArr) {
//                                        Transcript[] parentTranscripts = e.getParentTranscript();
//                                        allTranscripts.addAll(Arrays.asList(parentTranscripts));
//                                        for (Transcript ts : parentTranscripts) {
//                                            Integer counter = transcriptCounter.get(ts);
//                                            if (counter == null) {
//                                                counter = 0;
//                                            }
//                                            counter++;
//                                            transcriptCounter.put(ts, counter);
//                                        }
//                                    }
//
//                                    for (Transcript ts : allTranscripts) {
//                                        Integer counter = transcriptCounter.get(ts);
//                                        if (counter.equals(mathedExons.size())) {
//                                            // this is one of the transcripts that is shared by the set of exons..
//                                        }
//                                    }
//                                }
//                            }
//                            } else {
//                                for (Exon e : exonsForTranscript) {
//                                    int exonStart = e.getStart();
//                                    int exonStop = e.getEnd();
//                                    if (chrpos >= exonStart && chrpos <= exonStop) {
//                                        probeMapsToTheseExons.add(e);
//                                    }
//                                }
//                            }
//		    for (Map.Entry<String, Gene> geneset : availableGenes) {
//
//			Gene gene = geneset.getValue();
//			int start = gene.getStart();
//			int end = gene.getEnd();
//
//                        ArrayList<Transcript> selectedTranscripts = new ArrayList<Transcript>();
//                        
////			if (chrpos >= start && chrpos <= end) {
////			    // System.out.println(i + "\t" + gene.getName());
////                            HashMap<String, Transcript> transcripts = gene.getTranscripts();
////                            Set<Entry<String, Transcript>> transcriptset = transcripts.entrySet();
////                            
////                            
////                            // 
////			    if (ensemblGenes.length() > 0) {
////				ensemblGenes += "," + gene.getName();
////				strands += "," + gene.getStrand();
////				TStart += "," + gene.getStart();
////				TEnd += "," + gene.getEnd();
////			    } else {
////				ensemblGenes += "" + gene.getName();
////				strands += "" + gene.getStrand();
////				TStart += "" + gene.getStart();
////				TEnd += "" + gene.getEnd();
////			    }
////
////			    if (hgnc.length() > 0) {
////				hgnc += "," + gene.getAnnotation();
////			    } else {
////				hgnc += "" + gene.getAnnotation();
////			    }
////
////
////
////			    ctr++;
////			}
//                        
//		    }
//			if(ctr == 1){
//                    if (ensemblGenes.length() == 0) {
//                        out.writeln(probes[i] + "\t" + oldSeq.get(probes[i]) + "\t" + ChrAnnotation.parseByte(chr) + "\t" + originalChrPos + "\t-\t-\t-\t-\t-");
//                    } else {
//                        out.writeln(probes[i] + "\t" + oldSeq.get(probes[i]) + "\t" + ChrAnnotation.parseByte(chr) + "\t" + originalChrPos + "\t" + hgnc + "\t" + ensemblGenes + "\t" + strands + "\t" + TStart + "\t" + TEnd);
//                    }
//
////			}
////		    System.out.println(probes[i] + "\t" + ChrAnnotation.parseByte(chr) + "\t" + originalChrPos + "\t" + ctr + "\t" + hgnc + "\t" + ensemblGenes);
//
//                } else {
////			System.out.println(i + "\t" + ctr);
//                    out.writeln(probes[i] + "\t" + oldSeq.get(probes[i]) + "\t-1\t-1\t-1\t-\t-\t-\t-\t-");
//                }
//            } else {
//                out.writeln(probes[i] + "\t" + oldSeq.get(probes[i]) + "\t-1\t-1\t-1\t-\t-\t-\t-\t-");
////		    System.out.println(i + "\t" + ctr);
//            }
////	    }
//
//        }
//
//        out.close();
//            }
    private void convertToProbeTranslationFile(String probeTranslationFile, String infile, String outfile) throws IOException {

        System.out.println(probeTranslationFile);
        System.out.println(infile);
        TextFile tf = new TextFile(probeTranslationFile, TextFile.R);

        String[] header = tf.readLineElems(TextFile.tab);

        ArrayList<String> availableAnnotations = new ArrayList<String>();
        for (int h = 5; h < header.length; h++) {
            availableAnnotations.add(header[h]);
            System.out.println("Added annotation:\t" + header[h]);
        }

        HashMap<String, String> probeToAnnotation = new HashMap<String, String>();

        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            String strAnnot = Strings.concat(elems, Strings.tab, 5, elems.length);
            probeToAnnotation.put(elems[0], strAnnot);

            elems = tf.readLineElems(TextFile.tab);
        }

        tf.close();

        TextFile tf2 = new TextFile(infile, TextFile.R);
        TextFile tfOut = new TextFile(outfile, TextFile.W);

        String head = "Probe\tSeq\tChr\tPos\tHUGO";

        for (String s : availableAnnotations) {
            head += "\t" + s;
        }

        tfOut.writeln(head);
        elems = tf2.readLineElems(TextFile.tab);
        elems = tf2.readLineElems(TextFile.tab);
        while (elems != null) {

            String pb = null;
            if (elems.length < 5) {
                pb = Strings.concat(elems, Strings.tab, 0, 4);
                pb += "\t-";
            } else {
                pb = Strings.concat(elems, Strings.tab, 0, 5);
            }

            String annotation = probeToAnnotation.get(elems[0]);
            String outln = pb + "\t" + annotation;

            tfOut.writeln(outln);
            elems = tf2.readLineElems(TextFile.tab);
        }


        tfOut.close();
        tf2.close();
    }

    private void convertToPlatformAnnotationFiles(String infile, String outdir) throws IOException {
        TextFile tf = new TextFile(infile, TextFile.R);

//	HashMap<String, String> probeToPlatformMap = new HashMap<String, String>(); // maps platform + probeNr -> array address
        String[] header = tf.readLineElems(TextFile.tab); // header


        tf.close();


        for (int i = 5; i < header.length; i++) { // for each annotation
            tf.open();
            TextFile outfile = new TextFile(outdir + "/2012-03-07-" + header[i] + "Annotation.txt.gz", TextFile.W);

            String outheeader = "platform\tarrayaddress\tchr\tstart\tend\thugo";
            outfile.writeln(outheeader);
            String[] elems = tf.readLineElems(TextFile.tab);
            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {


                String probe = elems[0];
                String chr = elems[2];
                String pos = elems[3];
                String hugo = elems[4];

                String start = "-1";
                String end = "-1";
                if (!pos.equals("-1")) {
                    String[] posElems = pos.split(":");
                    String[] posElems2 = posElems[0].split("-");
                    String[] posElems3 = posElems[posElems.length - 1].split("-");

                    start = posElems2[0];
                    end = posElems3[posElems3.length - 1];
                }

                if (elems.length > i) {
                    String out = header[i] + "\t" + elems[i] + "\t" + hugo + "\t" + chr + "\t" + start + "\t" + end;
                    outfile.writeln(out);
                }


                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            outfile.close();
        }







    }
}
