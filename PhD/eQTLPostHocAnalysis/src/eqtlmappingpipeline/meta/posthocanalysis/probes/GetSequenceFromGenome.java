/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.probes;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class GetSequenceFromGenome {

    public static void main(String[] args) {
	GetSequenceFromGenome g = new GetSequenceFromGenome();
	try {
	    System.out.println(g.getPositionFromGenome(21, 43063890, 43063940));
	    System.out.println(g.getPositionFromGenome(21, 43068519, 43068569));
	    System.out.println(g.getPositionFromGenome(21, 43068601, 43068651));
	    System.out.println("CGCGGGATCCTTGTGCAGGGAAGAGCTGCCCTGGGCACCTGGCACCACAA");
	} catch (IOException ex) {
	    Logger.getLogger(GetSequenceFromGenome.class.getName()).log(Level.SEVERE, null, ex);
	}

    }

    /*
     * Platform probe transcript chr chrstart chrend strand matches mismatches gapQ gapT uniquehits 43063891	43063940 43068520	43068569 43068602	43068651
     */
    private String getPositionFromGenome(Integer chr, Integer left, Integer right) throws IOException {
	TextFile chromosomeSequenceFile = new TextFile("/Data/GenomeSequences/Ensembl54_HG18/dna/Homo_sapiens.NCBI36.54.dna.chromosome." + chr + ".fa.gz", TextFile.R);

	String line = chromosomeSequenceFile.readLine();
	int basesRead = 0;

	StringBuilder seq = new StringBuilder();
	while (line != null) {
	    if (line.startsWith(">")) {
		// fasta header line
	    } else {


		for (int i = 0; i < line.length(); i++) {

		    char c = line.charAt(i);
		    if (c == 'A' || c == 'T' || c == 'G' || c == 'C' || c == 'N') {
			if (basesRead >= left && basesRead < right) {
			    seq.append(c);
			} else if (basesRead >= left && basesRead >= right) {
			    break;
			}

		    } else {
			System.err.println("Could not parse char: " + c);
		    }
		    basesRead++;
		}
	    }

	    line = chromosomeSequenceFile.readLine();
	}
	chromosomeSequenceFile.close();
	return seq.toString();
    }
}
