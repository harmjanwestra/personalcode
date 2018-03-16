/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.zscoretables;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class zScoreTableSplitter {

    private GWASCatalog catalog;
    private ProbeTranslation probeTranslation;

    public void run(String gwasCatalog, String probetranslation, String zscoretable, String out) throws IOException, Exception {
	catalog = new GWASCatalog();
	catalog.read(gwasCatalog);

	probeTranslation = new ProbeTranslation();
	probeTranslation.load(probetranslation);

	if (!out.endsWith("/")) {
	    out += "/";
	}
	Gpio.createDir(out);


	HashSet<String> snpsInDataset = new HashSet<String>();
	TextFile zscoreFile = new TextFile(zscoretable, TextFile.R);
	int numSNPs = zscoreFile.countLines();
	int numProbes = zscoreFile.countCols(TextFile.tab) - 3;
	zscoreFile.close();
	zscoreFile.open();


	Float[][] zscores = new Float[numSNPs][numProbes];
	HashMap<String, Integer> snpNameToRow = new HashMap<String, Integer>();

	System.out.println("Parsing ZScoreFile\t" + zscoretable);

	// parse the header
	String[] columnnames = zscoreFile.readLineElemsReturnReference(TextFile.tab);
	String[] probenames = new String[numProbes];
	for (int i = 3; i < columnnames.length; i++) {
	    probenames[i - 3] = columnnames[i];
	}

// start parsing SNP content
	String[] zscoreelems = zscoreFile.readLineElemsReturnObjects(TextFile.tab);
	int lnctr = 0;
	ProgressBar pb = new ProgressBar(numSNPs);
	while (zscoreelems != null) {
	    String snpName = zscoreelems[0];
	    snpsInDataset.add(snpName);
	    snpNameToRow.put(snpName, lnctr);
	    String alleleAssessed = zscoreelems[2];
	    for (int i = 3; i < zscoreelems.length; i++) {
		Float zscore = null;
		try {
		    zscore = Float.parseFloat(zscoreelems[i]);
		} catch (NumberFormatException e) {
		}
		zscores[lnctr][i - 3] = zscore;
	    }

	    zscoreelems = null;
	    zscoreelems = zscoreFile.readLineElemsReturnObjects(TextFile.tab);

	    lnctr++;
	    pb.set(lnctr);
	}
	pb.close();

	System.out.println(snpsInDataset.size() + " SNPs have been meta-analysed");
	zscoreFile.close();

	GWASTrait[] traits = catalog.getTraits();
	for (GWASTrait t : traits) {

	    zscoreFile = new TextFile(zscoretable, TextFile.R);
	    GWASSNP[] snps = t.getSNPs();
	    HashSet<String> snpsForTrait = new HashSet<String>();

	    HashMap<String, String> alleleAssociatedWithTrait = new HashMap<String, String>();

	    for (GWASSNP s : snps) {
		String snpname = s.getName();

		if (snpsInDataset.contains(snpname)) {
		    snpsForTrait.add(snpname);
		    String riskAllele = s.getRiskAllele(t);
		    alleleAssociatedWithTrait.put(snpname, riskAllele);
		}
	    }



	    System.out.println(t.getName() + " - " + snpsForTrait.size());

	    if (snpsForTrait.size() > 0) {
		int numsnps = snpsForTrait.size();

		Float[][] data = new Float[numsnps][numProbes];
		String[] snpnames = new String[numsnps];

		for (int i = 0; i < snpnames.length; i++) {
		    String snp = snpnames[i];
		    Integer snpPos = snpNameToRow.get(snp);
		    data[i] = zscores[snpPos];
		}

		String traitname = t.getName();

		traitname = traitname.replaceAll("/", "_");
		traitname = traitname.replaceAll("\\\\", "_");
		traitname = traitname.replaceAll(" ", "_");

		TextFile outFile = new TextFile(out + numsnps + "-" + traitname + ".txt.gz", TextFile.W);
		String outheader = "-";
		for (int snp = 0; snp < numsnps; snp++) {
		    outheader += "\t" + snpnames[snp];
		}
		outFile.writeln(outheader);

		for (int i = 0; i < numProbes; i++) {
		    String line = null;
		    boolean probeHasNulls = false;
		    for (int j = 0; j < numsnps; j++) {
			if (data[j][i] == null) {
			    probeHasNulls = true;
			}
		    }

		    if (!probeHasNulls) {
			line = probenames[i];
			for (int j = 0; j < numsnps; j++) {
			    line += "\t" + data[j][i];
			}
			outFile.writeln(line);
		    }
		}

		outFile.close();
	    }


	}

    }
}
