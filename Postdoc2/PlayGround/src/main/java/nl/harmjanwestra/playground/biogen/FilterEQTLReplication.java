package nl.harmjanwestra.playground.biogen;

import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;

public class FilterEQTLReplication {
	
	public static void main(String[] args) {
		String input = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\brainseq\\allEqtls-rewrite.txt.gz";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\brainseq\\allEqtls-rewritegenes.txt.gz";
		String annotation = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\annotation\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz";
		
		FilterEQTLReplication f = new FilterEQTLReplication();
		try {
			f.rewriteGeneNames(input, annotation, output);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void filterFDR(String in, double threshold, String out) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfo = new TextFile(out, TextFile.W);
		tfo.writeln(tf.readLine());
		
		
		tf.close();
		tfo.close();
		
	}
	
	public void rewriteGeneNames(String in, String geneannot, String out) throws IOException {
		
		GTFAnnotation annot = new GTFAnnotation(geneannot);
		Collection<Gene> genes = annot.getGenes();
		HashMap<String, String> oldIdToNewID = new HashMap<String, String>();
		for (Gene gene : genes) {
			String geneName = gene.getName();
			String oldid = geneName.split("\\.")[0];
			oldIdToNewID.put(oldid, geneName);
		}
		
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfo = new TextFile(out, TextFile.W);
		String header = tf.readLine();
		tfo.writeln(header);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		int written = 0;
		while (elems != null) {
			String gene = elems[4];
			String newId = oldIdToNewID.get(gene);
			if (newId != null) {
				elems[4] = newId;
				tfo.writeln(Strings.concat(elems, Strings.tab));
				written++;
			}
			ctr++;
			if (ctr % 10000 == 0) {
				System.out.print("\r" + ctr + " lines parsed, " + written + " written.");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfo.close();
		System.out.println("Done");
		
	}
}
