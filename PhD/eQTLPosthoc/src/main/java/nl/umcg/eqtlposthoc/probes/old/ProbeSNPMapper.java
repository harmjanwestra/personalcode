/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.probes.old;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ProbeSNPMapper {

    ProbeTranslation p = null;
    HashMap<Byte, HashSet<Integer>> chromosomePositionsWithSNPs = new HashMap<Byte, HashSet<Integer>>();
    HashMap<Byte, HashMap<Integer, String>> chromosomeSNPPositions = new HashMap<Byte, HashMap<Integer, String>>();

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	ProbeSNPMapper s = new ProbeSNPMapper();
	try {


	    s.loadProbeAnnotation("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-02-22-ProbeAnnotation96PercIdentity-ExonJunctionMappingsFixed.txt");
//	    s.loadSNPAnnotation("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-27-SNPMappings-dbSNP130.txt.gz");
	    s.loadSNPAnnotation("/Data/GeneticalGenomicsDatasets/HapMap2r24-CEU/SNPMappings.txt");
	    s.mapSNPsToProbe("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-02-22-ProbeAnnotation96PercIdentity-ExonJunctionMappingsFixed.MappedAgainstHapMap2r24SNPs.txt");
//	    s.countTheNumberOfSNPsInTheVicinityOfProbes("", "");
	} catch (IOException e) {
	    e.printStackTrace();
	}
    }

    private void mapSNPsToProbe(String out) throws IOException {
	TextFile outfile = new TextFile(out, TextFile.W);
	String[] probes = p.getProbes();
	System.out.println("probe\tchr\tmidpos\tlengthofmapping\tnumsnps\tsnpswithin250kbofmidpos\tsnpsmapping");
	outfile.writeln("probe\tchr\tlengthofmapping\tnumsnps\tsnpswithin250kbofmidpos\tsnpsmapping");
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
	    outfile.writeln(i + "\t" + chr + "\t" + chrpos+ "\t"+ midpos + "\t" + probelengthsum + "\t" + snpcounter + "\t" + snpvicinitycounter + "\t" + allSNPNames);

	    if (probelengthsum > 50 || (probelengthsum > 0 && probelengthsum < 50)) {
		System.out.println(i + "\t" + chr + "\t" + chrpos+ "\t"+ midpos + "\t" + probelengthsum + "\t" + snpcounter + "\t" + snpvicinitycounter + "\t" + allSNPNames);
	    } else {
		System.out.println(i + "\t-\t-\t"+ midpos + "\t" + probelengthsum + "\t-\t-\t-");
	    }


	}
	outfile.close();
    }

    private void countTheNumberOfSNPsInTheVicinityOfProbes(String annotation, String dbsnp) throws IOException {
	throw new UnsupportedOperationException("Not yet implemented");
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

    private void loadProbeAnnotation(String probes) throws IOException {
	p = new ProbeTranslation();
	p.load(probes);
    }
}
