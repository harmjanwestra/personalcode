package nl.harmjanwestra.playground.conv;


import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class TriTyperSampleFilter {

	public static void main(String[] args) {

//        String gtdir = "D:\\geuvadis\\genotypes\\";
//        String gtoutdir = "D:\\geuvadis\\genotypes-EUR\\";
//        String samplefilter = "D:\\geuvadis\\GD660.GeneQuantCount-EUR-samples.txt";
//

		if (args.length < 3) {
			System.out.println("Usage: TTdirin TTdirout listofsamplestoinclude");

		} else {
			String gtdir = args[0];
			String gtoutdir = args[1];
			String samplefilter = args[2];

			TriTyperSampleFilter t = new TriTyperSampleFilter();
			try {
				t.run(gtdir, gtoutdir, samplefilter);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	public void run(String gtdir, String gtoutdir, String samplefilter) throws IOException {
		Gpio.createDir(gtoutdir);

		HashSet<String> incSamples = null;
		if (samplefilter != null) {
			TextFile tf = new TextFile(samplefilter, TextFile.R);
			ArrayList<String> list = tf.readAsArrayList();
			tf.close();
			incSamples = new HashSet<>();
			incSamples.addAll(list);
		}

		TextFile tf = new TextFile(gtdir + "Individuals.txt", TextFile.R);
		ArrayList<String> indlist = tf.readAsArrayList();
		tf.close();

		ArrayList<String> overlap = new ArrayList<>();
		TextFile ofi = new TextFile(gtoutdir + "Individuals.txt", TextFile.W);
		boolean[] inccol = new boolean[indlist.size()];
		for (int i = 0; i < indlist.size(); i++) {
			String ind = indlist.get(i);

			if (incSamples.contains(ind)) {
				overlap.add(ind);
				inccol[i] = true;
				ofi.writeln(ind);
			}
		}
		ofi.close();

		System.out.println("overlap: " + overlap.size());


		Gpio.copyFile(gtdir + "SNPMappings.txt.gz", gtoutdir + "SNPMappings.txt.gz");
		Gpio.copyFile(gtdir + "PhenotypeInformation.txt", gtoutdir + "PhenotypeInformation.txt");

		TextFile tfs = new TextFile(gtdir + "SNPs.txt.gz", TextFile.R);

		TextFile tfo = new TextFile(gtoutdir + "SNPs.txt.gz", TextFile.W);
		BinaryFile bfi = new BinaryFile(gtdir + "GenotypeMatrix.dat", BinaryFile.R);
		BinaryFile bfo = new BinaryFile(gtoutdir + "GenotypeMatrix.dat", BinaryFile.W);
		BinaryFile bfimp = null;
		BinaryFile bfomp = null;
		if (Gpio.exists(gtdir + "ImputedDosageMatrix.dat")) {
			bfimp = new BinaryFile(gtdir + "ImputedDosageMatrix.dat", BinaryFile.R);
			bfomp = new BinaryFile(gtoutdir + "ImputedDosageMatrix.dat", BinaryFile.R);
		}


		String ln = tfs.readLine();
		byte[] bufferA = new byte[indlist.size()];
		byte[] bufferB = new byte[indlist.size()];
		byte[] bufferC = new byte[indlist.size()];
		byte[] outbufferA = new byte[overlap.size()];
		byte[] outbufferB = new byte[overlap.size()];
		byte[] outbufferC = new byte[overlap.size()];

		int nrlinesparsed = 0;
		while (ln != null) {
			bfi.read(bufferA);
			bfi.read(bufferB);

			if (bfimp != null) {
				bfimp.read(bufferC);
			}

			int ctr = 0;
			for (int i = 0; i < indlist.size(); i++) {
				if (inccol[i]) {
					outbufferA[ctr] = bufferA[i];
					outbufferB[ctr] = bufferB[i];
					if (bfimp != null) {
						outbufferC[ctr] = bufferC[i];
					}
					ctr++;
				}
			}

			bfo.write(outbufferA);
			bfo.write(outbufferB);
			if (bfomp != null) {
				bfomp.write(outbufferC);
			}
			tfo.writeln(ln);
			ln = tfs.readLine();
			nrlinesparsed++;
			if (nrlinesparsed % 10000 == 0) {
				System.out.println(nrlinesparsed + " lines parsed");
			}
		}

		bfi.close();
		bfo.close();
		if (bfimp != null) {
			bfimp.close();
			bfomp.close();
		}
		tfo.close();
		tfs.close();


	}
}
