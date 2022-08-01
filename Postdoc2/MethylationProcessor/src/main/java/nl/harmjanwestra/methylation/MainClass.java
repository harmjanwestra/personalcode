package nl.harmjanwestra.methylation;

import umcg.genetica.math.matrix2.DoubleMatrixConverter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRandomAccessTranspose;

import java.io.IOException;

public class MainClass {
    public static void main(String[] args) {
//
//		String[] exclusion = new String[]{
//				"Z:\\projects\\2018-methylation\\code\\450K_DataProcessing\\ADDITIONAL_INFO\\ProbeFiltering\\freq5percent\\probeToFilter_450K_1000G_omni2.5.hg19.EUR_alleleFreq5percent_50bp_wInterroSite.txt",
//				"Z:\\projects\\2018-methylation\\code\\450K_DataProcessing\\ADDITIONAL_INFO\\ProbeFiltering\\ProbesBindingNonOptimal\\Source&BSProbesMappingMultipleTimesOrNotBothToBSandNormalGenome.txt"
//		};
//
//		String infolder = "Z:\\projects\\2018-methylation\\GPL13534_450k_OUT\\";
//		String probelistout = "Z:\\projects\\2018-methylation\\probelists\\GPL13534_450k";
//		try {
//			m.makeprobelistPerType(infolder, probelistout, exclusion);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

        try {
            if (args.length == 0) {
                printUsage();
            } else if (args[0].equals("makeprobelist")) {
                if (args.length < 4) {
                    System.out.println("Usage: makeprobelistfile infolder probeoutdir exclusionlist");
                    System.out.println("Exclusion list comma separated list of files");
                } else {
                    Legacy l = new Legacy();
                    l.makeprobelistPerType(args[1], args[2], args[3].split(";"));
                }
            } else if (args[0].equals("mergetables")) {
                if (args.length < 5) {
                    System.out.println("Usage: mergetables probelist1 probelist2 infolder outfolder");
                } else {
                    Legacy l = new Legacy();
                    l.runappendIndependentMatrices(args[1], args[2], args[3], args[4]);
                }
            } else if (args[0].equals("qqnorm")) {
                if (args.length < 3) {
                    System.out.println("Usage: qqnorm input output");
                } else {
                    QuantileNormalizeDiskBased450K q = new QuantileNormalizeDiskBased450K();
                    q.run(args[1], args[2]);
                }
            } else if (args[0].equals("logtransform")) {
                if (args.length < 3) {
                    System.out.println("Usage: logtransform input output");
                } else {
                    Log2TransformDiskBased450K q = new Log2TransformDiskBased450K();
                    q.run(args[1], args[2]);

                }
            } else if (args[0].equals("countnans")) {
                if (args.length < 3) {
                    System.out.println("Usage: countnans input output");
                } else {
                    DetermineNAs450K d = new DetermineNAs450K();
                    d.run(args[1], args[2]);
                }
            } else if (args[0].equals("calcmvals")) {

                if (args.length < 4) {
                    System.out.println("Usage: calcmvals inU inM out");
                } else {
                    CalculateMValuesDiskBased450K d = new CalculateMValuesDiskBased450K();
                    d.run(args[1], args[2], args[3]);
                }
            } else if (args[0].equals("calcbvals")) {

                if (args.length < 4) {
                    System.out.println("Usage: calcbvals inU inM out");
                } else {
                    CalculateBetaValuesDiskBased450K d = new CalculateBetaValuesDiskBased450K();
                    d.run(args[1], args[2], args[3]);
                }
            }else if (args[0].equals("mergebinarymatrices")) {
                Merge450KDiskBased d = new Merge450KDiskBased();
                if (args.length < 4) {
                    System.out.println("Usage: inT1 inT2 out [probelist]");
                } else {
                    String pblist = null;
                    if (args.length >= 5) {
                        pblist = args[4];
                    }
                    d.run(pblist, args[1], args[2], args[3]);
                }
            } else if (args[0].equals("correl")) {
                CorrelateDiskBased450K d = new CorrelateDiskBased450K();
                if (args.length < 4) {
                    System.out.println("Usage: correl in out nrrowsperblock");
                } else {
                    d.runblocks(args[1], args[2], Integer.parseInt(args[3]));
                }
            } else if (args[0].equals("detriangle")) {
                if (args.length < 3) {
                    System.out.println("Usage: input output");
                } else {
                    DeTriangularize d = new DeTriangularize();
                    d.run(args[1], args[2]);
                }
            } else if (args[0].equals("pca")) {
                if (args.length < 3) {
                    System.out.println("Usage: in nrpcs");
                } else {
                    PCA p = new PCA();
                    p.run(args[1], Integer.parseInt(args[2]));
                }
            } else if (args[0].equals("filter")) {
                if (args.length < 7) {
                    System.out.println("Usage: input output (rowset|null) (include|exclude) (colset|null) (include|exclude)");
                } else {
                    DatasetFilter filter = new DatasetFilter();
                    filter.run(args[1], args[2], args[3], args[4], args[5], args[6]);
                }
            } else if (args[0].equals("centerscale")) {
                CenterScale450K c = new CenterScale450K();
                if (args.length < 3) {
                    System.out.println("Usage; input output");
                } else {
                    c.run(args[1], args[2]);
                }
            } else if (args[0].equals("center")) {
                CenterScale450K c = new CenterScale450K();
                if (args.length < 3) {
                    System.out.println("Usage; input output");
                } else {
                    c.runCenterOnly(args[1], args[2]);
                }
            } else if (args[0].equals("removezerovariancerows")) {
                CenterScale450K c = new CenterScale450K();
                if (args.length < 3) {
                    System.out.println("Usage; input output");
                } else {
                    c.removeZeroVariance(args[1], args[2]);
                }
            } else if (args[0].equals("correctcovariates")) {
                CovariateCorrect450K c = new CovariateCorrect450K();
                if (args.length < 4) {
                    System.out.println("Usage; input covariates output");
                } else {
                    c.correct(args[1], args[2], args[3]);
                }
            } else if (args[0].equals("txttobin")) {
                if (args.length < 3) {
                    System.out.println("Usage: input output");
                } else {
                    DoubleMatrixConverter.TextToBinary(args[1], args[2]);
                }
            } else if (args[0].equals("bintotxt")) {
                if (args.length < 3) {
                    System.out.println("Usage: input output");
                } else {
                    DoubleMatrixConverter.BinaryToText(args[1], args[2]);
                }
            } else if (args[0].equals("transpose")) {
                if (args.length < 3) {
                    System.out.println("Usage: input output [nrrowstostoreinmemory:1000]");
                } else {
                    DoubleMatrixDatasetRandomAccessTranspose transpose = new DoubleMatrixDatasetRandomAccessTranspose();
                    int nrRowsToProcess = 1000;
                    if (args.length == 4) {
                        nrRowsToProcess = Integer.parseInt(args[3]);
                    }
                    transpose.transposeLargeMatrix(args[1], args[2], nrRowsToProcess);

                }
            } else if (args[0].equals("compare")) {
                if (args.length < 3) {
                    System.out.println("Usage: input1 input2");
                } else {
                    CompareBinaryMatrices c = new CompareBinaryMatrices();
                    c.run(args[1], args[2]);
                }
            }  else {
                printUsage();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private static void printUsage() {
        System.out.println("Usage:\n" +
                "txttobin\tConvert text to binary matrix\n" +
                "bintotxt\tConvert binary to text matrix\n\n" +
                "makeprobelist\tlist probes in output (probes expected on rows, samples on cols)\n" +
                "mergetables\tmerge 450K pipeline folder output\n" +
                "calcmvals\tcalculate m-values\n" +
                "calcbvals\tcalculate beta-values\n" +
                "qqnorm\tquantile normalize (expect samples on rows)\n" +
                "logtransform\tlog2transform (expect samples on rows)\n" +
                "correl\tcorrelations (expect samples on rows); outputs upper triangle correlation matrix over rows\n" +
                "countnans\tcount the number of nans in a matrix\n" +
                "mergebinarymatrices\tmerge binary matrices (expect samples on rows)\n" +
                "detriangle\tfill in the lower triangle of the matrix\n" +
                "pca\tcalculate eigenvectors\n" +
                "centerscale\tcenter and scale data (expect samples on rows)\n" +
                "center\tonly center (rows)\n" +
                "correctcovariates\tcorrect for covariates\n" +
                "filter\tfilter rows or columns\n" +
                "transpose\ttranspose binary matrix\n" +
                "removezerovariancerows\tRemove zero variance rows\n" +
                "compare\tCompare two binary matrices\n");

    }


}
