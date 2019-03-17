package nl.harmjanwestra.playground;


import umcg.genetica.containers.Pair;
import umcg.genetica.io.trityper.ProbeAnnotation;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class FDRBenchmark {
	
	
	enum FDRMETHOD {
		GW, SNPGENE, GENEALLSNPS, SNPALLGENES
	}
	
	public void run(String ttdir, String exp, String eneannotation) throws Exception {
		DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(exp);
		
		TriTyperGenotypeData d = new TriTyperGenotypeData(ttdir);
		
		// link samples
		String[] indsGt = d.getIndividuals();
		ArrayList<String> indsEx = ds.getColObjects();
		
		HashSet<String> indsExHash = new HashSet<>();
		indsExHash.addAll(indsEx);
		
		ArrayList<String> shared = new ArrayList<>();
		for (String ind : indsGt) {
			if (indsExHash.contains(ind)) {
				shared.add(ind);
			}
		}
		
		// index
		HashMap<String, Integer> indsExMap = ds.getHashCols();
		int[] expInddIndex = new int[shared.size()];
		int[] gtIndIndex = new int[shared.size()];
		for (int i = 0; i < shared.size(); i++) {
			String ind = shared.get(i);
			expInddIndex[i] = indsExMap.get(ind);
			gtIndIndex[i] = d.getIndividualId(ind);
		}
		
		// load gene annotation
		ProbeAnnotation pa = new ProbeAnnotation();
		pa.load(eneannotation);
		
		//
		String[] annotatedgenes = pa.getProbes();
		ArrayList<Integer> chrgenes = new ArrayList<>();
		
		for (int p = 0; p < annotatedgenes.length; p++) {
			if (pa.getChr()[p] == 22) {
				String name = pa.getProbes()[p];
				if (ds.getHashRows().containsKey(name)) {
					Integer id = ds.getHashRows().get(name);
					chrgenes.add(id);
				}
				
			}
		}
		// perform analysis
		ArrayList<Pair<Integer, Integer>> testpairsSNPs = new ArrayList<Pair<Integer, Integer>>();
		ArrayList<Pair<Integer, Integer>> testpairsGenes = new ArrayList<Pair<Integer, Integer>>();
		for (int s = 0; s < d.getSNPs().length; s++) {
			SNP obj = d.getSNPObject(s);
			if (obj.getChr() == 22) {
				for (int p = 0; p < chrgenes.size(); p++) {
					Integer id = chrgenes.get(p);
					if (Math.abs(pa.getChrStart()[id] - obj.getChrPos()) < 250000) {
						testpairsSNPs.add(new Pair<Integer, Integer>(s, id));
					}
				}
			}
		}
		
		SNPLoader loader = d.createSNPLoader();
		for (int s = 0; s < d.getSNPs().length; s++) {
			SNP obj = d.getSNPObject(s);
			if (obj.getChr() == 22) {
			
			}
		}
		
	}
}
