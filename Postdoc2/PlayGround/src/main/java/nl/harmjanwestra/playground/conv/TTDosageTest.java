package nl.harmjanwestra.playground.conv;

import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

import java.io.IOException;

public class TTDosageTest {

	public static void main(String[] args) {

		try {
			TriTyperGenotypeData d = new TriTyperGenotypeData(args[0]);
			SNPLoader loader = d.createSNPLoader();
			for (int s = 0; s < d.getSNPs().length; s++) {
				SNP obj = d.getSNPObject(s);
				loader.loadGenotypes(obj);
				loader.loadDosage(obj);

				if (obj.getMAF() > 0.25 && obj.getHWEP() > 0.0001) {
					for (int i = 0; i < d.getIndividuals().length; i++) {
						System.out.println(i + "\t" + obj.getGenotypes()[i] + "\t" + obj.getDosageValues()[i]);
					}
					loader.close();
					System.exit(-1);
				}

			}
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
