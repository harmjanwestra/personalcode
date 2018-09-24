package nl.harmjanwestra.playground.methylation;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

public class TransposeMatrix {
	
	public static void main(String[] args) {
		TransposeMatrix t = new TransposeMatrix();
		try {

//            t.run(args[0], args[1]);
			t.run("D:\\TMP\\methylation\\test.txt", "D:\\TMP\\methylation\\test-t.txt");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void run(String in, String out) throws Exception {

//        System.out.println(in);
//        System.out.println(out);
//        TextFile tf = new TextFile(in, TextFile.R);
//        TextFile tfo = new TextFile(in + "-tmp.txt", TextFile.W);
//        tfo.writeln(tf.readLine());
//        String ln = tf.readLine();
//        HashSet<String> rowsseen = new HashSet<String>();
//        int ctr = 0;
//        while (ln != null) {
//            String[] elems = Strings.subsplit(ln, Strings.tab, 0, 1);
//            if (!rowsseen.contains(elems[0])) {
//                tfo.writeln(ln);
//                rowsseen.add(elems[0]);
//            }
//            ctr++;
//            if (ctr % 100 == 0) {
//                System.out.print("\r"+ctr + " lines read");
//            }
//            ln = tf.readLine();
//        }
//        tf.close();
//        tfo.close();
//        System.out.println();
//        System.out.println("Moving " + in + "-tmp.txt.gz to " + in);
//        Gpio.moveFile(in + "-tmp.txt.gz", in);
		
		DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(in);
		DoubleMatrix2D matrix = ds.getMatrix();
		DenseDoubleAlgebra d = new DenseDoubleAlgebra();
		DoubleMatrix2D matout = d.transpose(matrix);
		
		DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<String, String>();
		dsout.setMatrix(matout);
		dsout.setColObjects(ds.getRowObjects());
		dsout.setRowObjects(ds.getColObjects());
		dsout.save(out);
	}
}
