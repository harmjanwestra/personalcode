package nl.harmjanwestra.playground.biogen.sqtl;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.containers.Pair;
import umcg.genetica.features.Exon;
import umcg.genetica.features.Gene;
import umcg.genetica.features.Transcript;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

public class AnnotateSQTLExons {

	public static void main(String[] args) {
//		String gtf = "D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\gencode.v32.primary_assembly.annotation.gtf.gz";
		String gtf = "D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz";

		String[] fin = new String[]{
				"D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\A3SS-merged-withqval.txt",
				"D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\A5SS-merged-withqval.txt",
				"D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\SE-merged-withqval.txt",
				"D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\RI-merged-withqval.txt",
				"D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\MXE-merged-withqval.txt"
		};

		String[] fout = new String[]{
				"D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\A3SS-merged-withqval-withexonnrs.txt",
				"D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\A5SS-merged-withqval-withexonnrs.txt",
				"D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\SE-merged-withqval-withexonnrs.txt",
				"D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\RI-merged-withqval-withexonnrs.txt",
				"D:\\Sync\\TMP\\sqtl\\qc\\12-spliceqtl\\MXE-merged-withqval-withexonnrs.txt"
		};

		AnnotateSQTLExons s = new AnnotateSQTLExons();
		try {
			s.run(fin, gtf, fout);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String[] filesIn, String gtfFile, String[] filesOut) throws IOException {
		GTFAnnotation gtf = new GTFAnnotation(gtfFile);
		System.out.println("Done reading GTF");


		for (int f = 0; f < filesIn.length; f++) {
			String fileIn = filesIn[f];
			System.out.println(fileIn);
			String fileOut = filesOut[f];
			TextFile tf = new TextFile(fileIn, TextFile.R);
			TextFile tfo = new TextFile(fileOut, TextFile.W);
			tfo.writeln(tf.readLine() + "\tStrand\tExonNumber");
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				String id = elems[0];
				String[] idElems = id.split("\\|");
				String gene = idElems[0];
				String[] posElems = idElems[1].split("_");
				ArrayList<Pair<Integer, Integer>> posPairs = new ArrayList<Pair<Integer, Integer>>();
				for (int q = 0; q < posElems.length; q += 2) {
					int p1 = Integer.parseInt(posElems[q]);
					int p2 = Integer.parseInt(posElems[q + 1]);
					posPairs.add(new Pair<>(p1, p2));
				}

				int stop = 1;
				ArrayList<ArrayList<Integer>> ranks = new ArrayList<>();
				ranks.add(new ArrayList<>());
				if (fileIn.contains("MXE")) {
					stop = 2;
					ranks.add(new ArrayList<>());
				}


				Gene geneobj = gtf.getStrToGene().get(gene);
				if (geneobj != null) {
					for (Transcript t : geneobj.getTranscripts()) {
						for (int pp = 0; pp < stop; pp++) {
							Pair<Integer, Integer> p = posPairs.get(pp);

							for (int er = 0; er < t.getExons().size(); er++) {
								Exon e = t.getExons().get(er);
								if (e.getStart() == p.getLeft() || e.getStop() == p.getRight()) {
									// matching exon
									Integer rank = t.getExonRank(e);
									System.out.println(pp + "\t" + p + "\t" + geneobj.getStrand() + "\t" + geneobj.getName() + "\t" + t.getName() + "\t" + e.getName() + "\t" + rank);

									ArrayList r = ranks.get(pp);
									r.add(rank);
								}
							}
						}
					}
				}

				String rankStr = "-";
				if (!ranks.get(0).isEmpty()) {
					rankStr = Strings.concat(ranks.get(0).toArray(new Integer[0]), Strings.comma);
				}

				if (fileIn.contains("MXE")) {
					if (!ranks.get(1).isEmpty()) {
						rankStr += ";" + Strings.concat(ranks.get(1).toArray(new Integer[0]), Strings.comma);
					} else {
						rankStr += ";-";
					}
				}

				tfo.writeln(Strings.concat(elems, Strings.tab) + "\t" + geneobj.getStrand() + "\t" + rankStr);

//			System.exit(0);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			tfo.close();
		}


	}
}
