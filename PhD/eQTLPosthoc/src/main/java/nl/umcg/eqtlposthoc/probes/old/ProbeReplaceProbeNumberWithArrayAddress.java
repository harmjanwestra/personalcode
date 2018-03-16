/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.probes.old;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ProbeReplaceProbeNumberWithArrayAddress {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here
	try {
	    ProbeReplaceProbeNumberWithArrayAddress pb = new ProbeReplaceProbeNumberWithArrayAddress();
	    pb.run("/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-01-11-MetaAnalysis/cis/2012-02-22-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/eQTLProbesFDR0.05.txt", "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-14-ProbeTranslationTable+H8HT12Conversion.txt");
	} catch (IOException e) {
	    e.printStackTrace();
	}
    }

    public void run(String in, String pt) throws IOException {

	TextFile prt = new TextFile(pt, TextFile.R);
	String match = "H8v2ConvToHT12";
	String[] headerelems = prt.readLineElems(TextFile.tab);

	int platformcol = 0;
	for (int i = 0; i < headerelems.length; i++) {
	    if (headerelems[i].contains(match)) {
		platformcol = i;
		break;
	    }
	}

	System.out.println("Platform col " + platformcol);

	HashMap<String, String> probeToArrayAddress = (HashMap<String, String>) prt.readAsHashMap(0, platformcol);

	System.out.println("probes loaded " + probeToArrayAddress.size());
	prt.close();

	TextFile eqtlin = new TextFile(in, TextFile.R);

	TextFile eqtlout = new TextFile(in + "AnnotatedFor" + match + ".txt", TextFile.W);

	String[] elems = eqtlin.readLineElems(TextFile.tab);
	String output = Strings.concat(elems, Strings.tab);
	eqtlout.writeln(output);

	elems = eqtlin.readLineElems(TextFile.tab);
	int inln = 0;
	int outln = 0;
	while (elems != null) {

	    String probe = elems[4];

	    String newprobe = probeToArrayAddress.get(probe);


	    inln++;
	    if (!newprobe.equals("-")) {
		elems[4] = newprobe;
		output = Strings.concat(elems, Strings.tab);
		eqtlout.writeln(output);
		outln++;
	    }



	    elems = eqtlin.readLineElems(TextFile.tab);
	}

	eqtlin.close();
	eqtlout.close();

	System.out.println(inln+" - "+ outln);

    }
}
