//package nl.harmjanwestra.playground.biogen.freeze2;
//
//import com.jujutsu.tsne.TSne;
//import umcg.genetica.io.text.TextFile;
//import umcg.genetica.math.matrix2.DoubleMatrixDataset;
//import umcg.genetica.math.stats.Descriptives;
//
//public class RNASeqOutliers {
//	public static void main(String[] args) {
//
////		String filesfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\evfiles.txt";
////		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\";
//		try {
////			TextFile tf = new TextFile(filesfile, TextFile.R);
////			String[] elems = tf.readLineElems(TextFile.tab);
////			while (elems != null) {
////				String ds = elems[0];
////				String dsname = elems[1];
////				RNASeqOutliers s = new RNASeqOutliers();
////				s.run(ds,
////						out + "RNA-zlt3-outliers-" + dsname + ".txt",
////						out + "RNA-zlt3-nonoutliers-" + dsname + ".txt");
////				elems = tf.readLineElems(TextFile.tab);
////			}
////			tf.close();
////
//			RNASeqOutliers s = new RNASeqOutliers();
////			s.run("D:\\Freeze2\\2019-04-11-Freeze2.TMM.SampleFilter.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors-pc1to4.txt",
////					"D:\\Freeze2\\2019-04-11-Freeze2.TMM.SampleFilter.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors-pc1to4-outliers.txt",
////					"D:\\Freeze2\\2019-04-11-Freeze2.TMM.SampleFilter.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors-pc1to4-nonoutliers.txt");
//
////			s.run("D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\2-samplefilter\\2019-04-11-Freeze2.TMM.SampleFilter.SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors-pc1to4.txt",
////					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\2-samplefilter\\2019-04-11-Freeze2.TMM.SampleFilter.SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors-pc1to4-outliers.txt",
////					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\2-samplefilter\\2019-04-11-Freeze2.TMM.SampleFilter.SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors-pc1to4-nonoutliers.txt");
//
//			s.run("D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\4-samplefilterrun2\\2019-04-11-Freeze2.TMM.SampleFilter.SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors-pc1to4.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\4-samplefilterrun2\\2019-04-11-Freeze2.TMM.SampleFilter.SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors-pc1to4-outliers.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\rnaqc\\2019-04-11\\4-samplefilterrun2\\2019-04-11-Freeze2.TMM.SampleFilter.SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors-pc1to4-nonoutliers.txt");
//
//		} catch (
//				Exception e) {
//			e.printStackTrace();
//		}
//
//	}
//
//	public void run(String file, String out, String out2) throws Exception {
//		DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(file);
//		double[] data = ds.viewCol(0).toArray();
//		double mean = Descriptives.mean(data);
//		double sd = Math.sqrt(Descriptives.variance(data));
//
//		TextFile tf = new TextFile(out, TextFile.W);
//		TextFile tf2 = new TextFile(out2, TextFile.W);
//		for (int i = 0; i < data.length; i++) {
//			double z = Math.abs((data[i] - mean) / sd);
//			if (z > 3) {
//				tf.writeln(ds.getRowObjects().get(i));
//			} else {
//				tf2.writeln(ds.getRowObjects().get(i));
//			}
//		}
//		tf.close();
//		tf2.close();
//
//	}
//}
