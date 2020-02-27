package nl.harmjanwestra.playground.biogen;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

public class ENSGToHGNC {

	public static void main(String[] args) {

		String input = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\annotation\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\annotation\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz";

		try {
			GTFAnnotation annot = new GTFAnnotation(input);
			TextFile tf = new TextFile("", TextFile.W);

			tf.writeln("ENSG\tHGNC");
			for (Gene g : annot.getGenes()) {
				tf.writeln(g.getName() + "\t" + g.getName());

			}

			tf.close();

		} catch (IOException e) {
			e.printStackTrace();
		}


	}


}
