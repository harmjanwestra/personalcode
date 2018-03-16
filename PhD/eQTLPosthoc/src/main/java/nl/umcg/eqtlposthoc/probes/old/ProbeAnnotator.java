/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.probes.old;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;
import umcg.genetica.containers.Chromosome;
import umcg.genetica.containers.Exon;
import umcg.genetica.containers.Gene;
import umcg.genetica.containers.Transcript;
import umcg.genetica.ensembl.Features;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harmjan
 */
public class ProbeAnnotator {

    public static void main(String[] args) {
        // TODO code application logic here

        try {

            ProbeAnnotator p = new ProbeAnnotator();
            p.run();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run() throws IOException {
        Features f = new Features();
        f.loadAnnotation("/Data/GenomeSequences/Ensembl54_HG18/structures/structures_b54.txt");

        ProbeTranslation p = new ProbeTranslation();
        p.load("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-02-22-ProbeAnnotation96PercIdentity-ExonJunctionMappingsFixed.txt");

        ChrAnnotation chran = new ChrAnnotation();

        TextFile out = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-02-22-ProbeAnnotation96PercIdentity-ExonJunctionMappingsFixed.EnsemblAnnotated.txt", TextFile.W);
        out.writeln("Probe\tChr\tPos\tCount\tHGNC\tEnsembls");

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

        for (int i = 0; i < numProbes; i++) {
//	    if (presentOnHT12v3[i]) {

            byte chr = p.getProbeChr(i);


            int ctr = 0;
            if (chr > 0 && chr < 25) {
                int chrpos = p.getProbeChrPos(i);
                String actualMappingPosition = p.getActualMappingPosition(i);
                if (chrpos > 0) {
                    Chromosome c = f.getChromosomeHash().get(ChrAnnotation.parseByte(chr));
                    if (c == null) {
                        System.out.println("Chromsome " + chr + " not found for probe " + i);
                    }
                    HashMap<String, Gene> genes = c.getGenesHash();
                    Set<Entry<String, Gene>> availableGenes = genes.entrySet();
                    String ensemblGenes = "";
                    String hgnc = "";

                    HashSet<Gene> genesMappingToRegion = new HashSet<Gene>();
                    for (Entry<String, Gene> geneset : availableGenes) {

                        Gene gene = geneset.getValue();
                        int start = gene.getStart();
                        int end = gene.getEnd();

                        if (chrpos >= start && chrpos <= end) {
                            genesMappingToRegion.add(gene);
                            ctr++;
                        }
                    }

//                    if (genesMappingToRegion.size() > 1) {
                    // see whether there are exons to which this probe is specific.
                    Gene[] genesArr = genesMappingToRegion.toArray(new Gene[0]);
                    HashSet<Exon> probeMapsToTheseExons = new HashSet<Exon>();
                    for (Gene g : genesArr) {
                        HashMap<String, Transcript> transcripts = g.getTranscripts();
                        Set<Entry<String, Transcript>> transcriptset = transcripts.entrySet();
                        for (Entry<String, Transcript> transcriptEntry : transcriptset) {
                            Transcript t = transcriptEntry.getValue();
                            Exon[] exonsForTranscript = t.getExonsRanked();
                            // first see whether the probe actually maps in an exon, and not accross..
                            for (Exon e : exonsForTranscript) {
                                int exonStart = e.getStart();
                                int exonStop = e.getEnd();
                                if (chrpos >= exonStart && chrpos <= exonStop) {
                                    probeMapsToTheseExons.add(e);
                                }
                            }
                        }
                        // now check whether the transcript mapped into an exon...
                    }

//                    } else {
//                        // System.out.println(i + "\t" + gene.getName());
//                        Gene gene = genesMappingToRegion.toArray(new Gene[0])[0];
//                        if (ensemblGenes.length() > 0) {
//                            ensemblGenes += "," + gene.getName();
//                        } else {
//                            ensemblGenes += "" + gene.getName();
//                        }
//
//                        if (hgnc.length() > 0) {
//                            hgnc += "," + gene.getAnnotation();
//                        } else {
//                            hgnc += "" + gene.getAnnotation();
//                        }
//
//                    }

//			if(ctr == 1){
                    String originalChrPos = p.getActualMappingPosition(i);
                    out.writeln(i + "\t" + ChrAnnotation.parseByte(chr) + "\t" + originalChrPos + "\t" + ctr + "\t" + hgnc + "\t" + ensemblGenes);
//			}
                    System.out.println(i + "\t" + ChrAnnotation.parseByte(chr) + "\t" + originalChrPos + "\t" + ctr + "\t" + hgnc + "\t" + ensemblGenes);

                } else {
//			System.out.println(i + "\t" + ctr);
                }
            } else {
//		    System.out.println(i + "\t" + ctr);
            }
//	    }

        }

        out.close();
    }
}
