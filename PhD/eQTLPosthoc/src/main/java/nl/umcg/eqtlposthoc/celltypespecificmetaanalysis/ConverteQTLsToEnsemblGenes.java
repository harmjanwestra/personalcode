/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import umcg.genetica.containers.Exon;
import umcg.genetica.containers.Gene;
import umcg.genetica.containers.Transcript;
import umcg.genetica.ensembl.Features;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ConverteQTLsToEnsemblGenes {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        String ensemblAnnotation = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/All.0.96%Identity-Merged-PerProbe-UniqueMappings-Ensembl70.txt";
        String eqtlfile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysis/Vector-7330546-2.txt";
        String ensemblfeatures = "/Volumes/iSnackHD/Data/GenomeSequences/Ensembl70_HG19/structures/structures_b70.txt";
        String probetranslationfile = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-ProbeAnnotationFile.txt";

        try {

            Features featureset = new Features();
            featureset.loadAnnotation(ensemblfeatures);

            ProbeTranslation pbt = new ProbeTranslation();
            HashMap<String, String> ht12v3ToProbe = new HashMap<String, String>();
            ht12v3ToProbe = pbt.getProbeTranslation(probetranslationfile, "HT12v3.txt", "Probe");

            HashMap<String, String> probeToHT12v3 = new HashMap<String, String>();
            probeToHT12v3 = pbt.getProbeTranslation(probetranslationfile, "Probe", "HT12v3.txt");

            // probe annotation..
            HashMap<String, String> probesMappingToEns = new HashMap<String, String>();
            HashMap<String, String> probeToChrStart = new HashMap<String, String>();
            HashMap<String, String> probeToChrStop = new HashMap<String, String>();
            HashMap<String, String> probeToChr = new HashMap<String, String>();




            TextFile annot = new TextFile(ensemblAnnotation, TextFile.R);
            TextFile tfOut = new TextFile(ensemblAnnotation + "-FullEnsemblAnnotation.txt", TextFile.W);

            String[] annotelems = annot.readLineElems(TextFile.tab);
            tfOut.writeln("MetaProbeId\tChr\tChrStart\tChrEnd\tHUGO\tHT12v3\tEnsemblGenes\tEnsemblTranscripts\tEnsemblExons");

            annot.readLineElems(TextFile.tab);
            while (annotelems != null) {
                String metaprobe = annotelems[0];
                String hugo = annotelems[7];
                String mapseq = annotelems[2];
                String chr = annotelems[3];
                String chrStaStr = annotelems[4];
                String chrStoStr = annotelems[5];
                String ens = annotelems[6];

                probeToChrStart.put(metaprobe, chrStaStr);
                probeToChrStop.put(metaprobe, chrStoStr);
                probeToChr.put(metaprobe, chr);

                if (mapseq.startsWith("ENSE")) {
                    probesMappingToEns.put(metaprobe, mapseq);
                    ens = mapseq;
                } else {
                    probesMappingToEns.put(metaprobe, ens);
                }



                String defaultOutput = metaprobe + "\t" + chr + "\t" + chrStaStr + "\t" + chrStoStr + "\t" + hugo;
                if (metaprobe == null || ens == null || ens.equals("-")) {
                    // meeh

                    tfOut.writeln(defaultOutput + "\t" + probeToHT12v3.get(metaprobe)
                            + "\t-"
                            + "\t-"
                            + "\t-");

                } else {

                    // handle exon mappings first
                    if (ens.startsWith("ENSE")) {
                        // we have our answer.. get the transcript for this set of exons..
                        String[] allExons = ens.split("\\+");
                        HashSet<String> allexonNames = new HashSet<String>();
                        allexonNames.addAll(Arrays.asList(allExons));

                        HashSet<Transcript> transcripts = new HashSet<Transcript>();
                        for (String exon : allExons) {
                            Exon x = featureset.getExonHash().get(exon);
                            if (x == null) {
                                System.err.println("ERROR: exon not found: " + exon);
                            } else {
                                Transcript[] parents = x.getParentTranscript();
                                if (parents == null) {
                                    System.err.println("ERROR: exon " + exon + " has no parents?");
                                } else {
                                    transcripts.addAll(Arrays.asList(parents));
                                }

                            }

                        }

                        HashSet<String> associatedTranscripts = new HashSet<String>();
                        HashSet<String> associatedGenes = new HashSet<String>();
                        for (Transcript t : transcripts) {
                            int exonctr = 0;
                            Exon[] associatedExons = t.getExonsRanked();
                            for (Exon x : associatedExons) {
                                if (allexonNames.contains(x.getName())) {
                                    exonctr++;
                                }
                            }
                            if (exonctr == allexonNames.size()) {
                                associatedTranscripts.add(t.getName());
                                associatedGenes.add(t.getParentGene().getName());
                            }
                        }



                        tfOut.writeln(defaultOutput + "\t" + probeToHT12v3.get(metaprobe)
                                + "\t" + Strings.concat(associatedGenes.toArray(new String[0]), Strings.semicolon)
                                + "\t" + Strings.concat(associatedTranscripts.toArray(new String[0]), Strings.semicolon)
                                + "\t" + Strings.concat(allexonNames.toArray(new String[0]), Strings.semicolon));

                    } else {
                        // probe maps to a (set of) gene(s)..

                        String[] enselems = ens.split(";");
                        HashSet<String> matchingTranscripts = new HashSet<String>();
                        HashSet<String> matchingExons = new HashSet<String>();
                        for (String gene : enselems) {
                            Gene g = featureset.getGeneHash().get(gene);
                            HashMap<String, Transcript> transcripts = g.getTranscripts();
                            Set<Map.Entry<String, Transcript>> set = transcripts.entrySet();

                            if (gene.equals("ENSG00000130402")) {
                                System.out.println("Found it: " + set.size());
                            }


                            HashSet<Transcript> transcriptList = new HashSet<Transcript>();
                            for (Map.Entry<String, Transcript> e : set) {
                                transcriptList.add(e.getValue());
                            }

                            Integer chrSta = Integer.parseInt(chrStaStr);
                            Integer chrSto = Integer.parseInt(chrStoStr);

                            for (Transcript t : transcriptList) {
                                Exon[] allExons = t.getExonsRanked();

                                for (Exon e : allExons) {
                                    if (chrSta >= e.getStart() && chrSto <= e.getEnd()) {
                                        matchingExons.add(e.getName());
                                        matchingTranscripts.add(t.getName());
                                    }
                                }
                            }
                        }

                        String transStr = Strings.concat(matchingTranscripts.toArray(new String[0]), Strings.semicolon);
                        String exonStr = Strings.concat(matchingExons.toArray(new String[0]), Strings.semicolon);

                        if (transStr.length() == 0) {
                            transStr = "-";
                        }
                        if (exonStr.length() == 0) {
                            exonStr = "-";
                        }

                        tfOut.writeln(defaultOutput + "\t" + probeToHT12v3.get(metaprobe)
                                + "\t" + Strings.concat(enselems, Strings.semicolon)
                                + "\t" + transStr
                                + "\t" + exonStr);

                    }
                }



                annotelems = annot.readLineElems(TextFile.tab);
            }
            annot.close();
            tfOut.close();

            TextFile output = new TextFile(eqtlfile + "-WithEnsemblAnnotation.txt", TextFile.W);
            TextFile tf = new TextFile(eqtlfile, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            output.writeln(Strings.concat(elems, Strings.tab) + "\tMetaProbeId\tEnsemblGenes\tEnsemblTranscripts\tEnsemblExons");
            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {

                String eqtl = elems[0];
                String Z = elems[1];
                String ht12probe = eqtl.split("-")[1];
                String metaprobe = ht12v3ToProbe.get(ht12probe);


                String ens = probesMappingToEns.get(metaprobe);
                if (metaprobe == null || ens == null || ens.equals("-")) {
                    // meeh
                    output.writeln(eqtl + "\t" + Z + "\t" + metaprobe
                            + "\t-"
                            + "\t-"
                            + "\t-");

                } else {

                    // handle exon mappings first
                    if (ens.startsWith("ENSE")) {
                        // we have our answer.. get the transcript for this set of exons..
                        String[] allExons = ens.split("\\+");
                        HashSet<String> allexonNames = new HashSet<String>();
                        allexonNames.addAll(Arrays.asList(allExons));

                        HashSet<Transcript> transcripts = new HashSet<Transcript>();
                        for (String exon : allExons) {
                            Exon x = featureset.getExonHash().get(exon);
                            if (x == null) {
                                System.err.println("ERROR: exon not found: " + exon);
                            } else {
                                Transcript[] parents = x.getParentTranscript();
                                if (parents == null) {
                                    System.err.println("ERROR: exon " + exon + " has no parents?");
                                } else {
                                    transcripts.addAll(Arrays.asList(parents));
                                }

                            }

                        }

                        HashSet<String> associatedTranscripts = new HashSet<String>();
                        HashSet<String> associatedGenes = new HashSet<String>();
                        for (Transcript t : transcripts) {
                            int exonctr = 0;
                            Exon[] associatedExons = t.getExonsRanked();
                            for (Exon x : associatedExons) {
                                if (allexonNames.contains(x.getName())) {
                                    exonctr++;
                                }
                            }
                            if (exonctr == allexonNames.size()) {
                                associatedTranscripts.add(t.getName());
                                associatedGenes.add(t.getParentGene().getName());
                            }
                        }



                        output.writeln(eqtl + "\t" + Z + "\t" + metaprobe
                                + "\t" + Strings.concat(associatedGenes.toArray(new String[0]), Strings.semicolon)
                                + "\t" + Strings.concat(associatedTranscripts.toArray(new String[0]), Strings.semicolon)
                                + "\t" + Strings.concat(allexonNames.toArray(new String[0]), Strings.semicolon));

                    } else {
                        // probe maps to a (set of) gene(s)..

                        String[] enselems = ens.split(";");
                        HashSet<String> matchingTranscripts = new HashSet<String>();
                        HashSet<String> matchingExons = new HashSet<String>();
                        for (String gene : enselems) {
                            Gene g = featureset.getGeneHash().get(gene);
                            HashMap<String, Transcript> transcripts = g.getTranscripts();
                            Set<Map.Entry<String, Transcript>> set = transcripts.entrySet();

                            if (gene.equals("ENSG00000130402")) {
                                System.out.println("Found it: " + set.size());
                            }


                            HashSet<Transcript> transcriptList = new HashSet<Transcript>();
                            for (Map.Entry<String, Transcript> e : set) {
                                transcriptList.add(e.getValue());
                            }

                            Integer chrSta = Integer.parseInt(probeToChrStart.get(metaprobe));
                            Integer chrSto = Integer.parseInt(probeToChrStop.get(metaprobe));

                            for (Transcript t : transcriptList) {
                                Exon[] allExons = t.getExonsRanked();

                                for (Exon e : allExons) {
                                    if (chrSta >= e.getStart() && chrSto <= e.getEnd()) {
                                        matchingExons.add(e.getName());
                                        matchingTranscripts.add(t.getName());
                                    }
                                }
                            }
                        }

                        String transStr = Strings.concat(matchingTranscripts.toArray(new String[0]), Strings.semicolon);
                        String exonStr = Strings.concat(matchingExons.toArray(new String[0]), Strings.semicolon);

                        if (transStr.length() == 0) {
                            transStr = "-";
                        }
                        if (exonStr.length() == 0) {
                            exonStr = "-";
                        }

                        output.writeln(eqtl + "\t" + Z + "\t" + metaprobe
                                + "\t" + Strings.concat(enselems, Strings.semicolon)
                                + "\t" + transStr
                                + "\t" + exonStr);

                    }
                }


                elems = tf.readLineElems(TextFile.tab);
            }
            output.close();
            tf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
