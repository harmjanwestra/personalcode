package nl.harmjanwestra.playground.cis;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.AssociationPanel;
import nl.harmjanwestra.utilities.graphics.panels.GenePanel;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class EQTLRegionalPlot {
	
	
	public static void main(String[] args) {
	
	}
	
	public void run(String gene, String eqtlfile, String ldfile, String geneannotation, String output) throws IOException, DocumentException {
		
		GTFAnnotation gtf = new GTFAnnotation(geneannotation);
		Gene geneObj = gtf.getStrToGene().get(gene);
		
		int midpoint = (geneObj.getStart() + geneObj.getStop()) / 2;
		Feature region = new Feature(geneObj.getChromosome(), midpoint - 1000000, midpoint + 1000000);
		
		
		ArrayList<Gene> genes = new ArrayList<Gene>();
		for (Gene g : gtf.getGenes()) {
			if (g.overlaps(region)) {
				genes.add(g);
			}
		}
		
		Grid g = new Grid(500, 500, 2, 2, 10, 10);
		
		AssociationPanel assocP = new AssociationPanel(1, 1);
		ArrayList<Triple<String, Integer, Double>> pvals = readEQTLFile(eqtlfile, gene);
		
		
		double[] ldvals = getLDData(ldfile, pvals);
		ArrayList<Pair<Integer, Double>> pvalssortedorwhatever = new ArrayList<Pair<Integer, Double>>();
		for (Triple<String, Integer, Double> q : pvals) {
			pvalssortedorwhatever.add(new Pair<>(q.getMiddle(), q.getRight()));
		}
		
		assocP.setDataSingleDs(region, null, pvalssortedorwhatever, "eQTLgen");
		assocP.setLDData(ldvals);
		g.addPanel(assocP);
		GenePanel gp = new GenePanel(1, 1);
		
		gp.setData(region, genes);
		g.addPanel(gp);
		
		g.draw(output);
		
	}
	
	private double[] getLDData(String ldfile, ArrayList<Triple<String, Integer, Double>> pvals) throws IOException {
		
		// get top position
		String maxsnp = null;
		int maxsnpindex = 0;
		Double maxp = 0d;
		int maxindex = 0;
		int ctr = 0;
		HashMap<String, Integer> snpindex = new HashMap<>();
		for (int q = 0; q < pvals.size(); q++) {
			Triple<String, Integer, Double> p = pvals.get(q);
			snpindex.put(p.getLeft(), ctr);
			
			if (p.getRight() > maxp) {
				maxp = p.getRight();
				maxindex = q;
				maxsnp = p.getLeft();
			}
			ctr++;
		}
		
		double[] ld = new double[pvals.size()];
		
		TextFile tf = new TextFile(ldfile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp1 = elems[0];
			String snp2 = elems[1];
			
			if (snp1.equals(maxsnp) || snp2.equals(maxsnp)) {
				Integer snpid = null;
				if (snp1.equals(maxsnp)) {
					snpid = snpindex.get(snp2);
				} else {
					snpid = snpindex.get(snp1);
				}
				double rsq = Double.parseDouble(elems[1]);
				ld[snpid] = rsq;
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		return ld;
	}
	
	private ArrayList<Triple<String, Integer, Double>> readEQTLFile(String eqtlfile, String gene) throws IOException {
		
		ArrayList<Triple<String, Integer, Double>> output = new ArrayList<>();
		
		TextFile tf = new TextFile(eqtlfile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			
			if (elems[4].equals(gene)) {
				Double pval = Double.parseDouble(elems[0]);
				
				String snp = elems[1];
				Integer pos = Integer.parseInt(elems[2]);
				
				Triple<String, Integer, Double> p = new Triple<>(snp, pos, -Math.log10(pval));
				output.add(p);
			}
			elems = tf.readLineElems(TextFile.tab);
			
		}
		tf.close();
		
		return output;
	}
}
