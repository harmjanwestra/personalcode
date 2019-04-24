package nl.harmjanwestra.playground.biogen.covariates;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.util.Arrays;
import java.util.HashSet;

public class CovariateFilter {


	public static void main(String[] args) {

		String[] set = new String[]{
				"INTRONIC_BASES",
				"PCT_INTRONIC_BASES",
				"PCT_MRNA_BASES",
				"PCT_CODING_BASES",
				"PCT_USABLE_BASES",
				"PCT_INTERGENIC_BASES",
				"INTERGENIC_BASES",
				"MEDIAN_3PRIME_BIAS",
				"PCT_UTR_BASES",
				"PCT_RIBOSOMAL_BASES",
				"MEDIAN_5PRIME_BIAS",
		};
		HashSet<String> s = new HashSet<String>();
		s.addAll(Arrays.asList(set));

		String file = "D:\\Work\\Freeze2\\run2\\2019-04-11-Freeze2.TMM.Covariates-Numeric.txt";
		String fileout = "D:\\Work\\Freeze2\\run2\\2019-04-11-Freeze2.TMM.Covariates-Numeric-Top10Covariates.txt";
		DoubleMatrixDataset<String, String> ds = null;
		try {
			ds = DoubleMatrixDataset.loadSubsetOfTextDoubleData(file, '\t', null, s);
			ds.save(fileout);
		} catch (Exception e) {
			e.printStackTrace();
		}


	}

}
