package nl.harmjanwestra.playground.cis.ld;


import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.IOException;
import java.util.ArrayList;

public class BinaryLDMatrixConverter {
	
	public static void main(String[] args) {
		
		String[] genes = new String[]{
				"ENSG00000105499",
				"ENSG00000134184",
				"ENSG00000197146",
				"ENSG00000213366"
		};
		
		BinaryLDMatrixConverter cv = new BinaryLDMatrixConverter();
		for (String g : genes) {
			try {
				cv.run("D:\\ld\\geneswithclearhits\\" + g + "-ldmatrix.bin", "D:\\ld\\geneswithclearhits\\" + g + "-ldmatrix.txt.gz");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		
	}
	
	public void run(String ldmatrixin, String ldmatrixout) throws IOException {
		
		
		System.out.println(ldmatrixin);
		BinaryFile bf = new BinaryFile(ldmatrixin, BinaryFile.R);
		
		int nrsnps = bf.readInt();
		ArrayList<String> snps = new ArrayList<String>();
		for (int i = 0; i < nrsnps; i++) {
			snps.add(bf.readString());
		}
		
		/*
		for (int i = 0; i < snpList.size(); i++) {
						for (int j = i + 1; j < snpList.size(); j++) {
							float f = (float) ldmatrix.getMatrix().getQuick(i, j);
							bf.writeFloat(f);
						}
					}
		 */
		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(nrsnps, nrsnps);
		for (int i = 0; i < nrsnps; i++) {
			ds.getMatrix().setQuick(i, i, 1d);
			for (int j = i + 1; j < nrsnps; j++) {
				float f = bf.readFloat();
				ds.getMatrix().setQuick(i, j, f);
				ds.getMatrix().setQuick(j, i, f);
			}
		}
		bf.close();
		
		ds.save(ldmatrixout);
		
		
	}
}