/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.celltypespecificmetaanalysis;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class CellTypeSpecificVectorIsolator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

//            String matrixIn = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/RotterdamStudy/MarjoleinVisit3/CellTypeSpecificityTestOutputMT2/CellTypeSpecificityMatrix.binary";
//            String query = "rs1247635-1570553";
//            String outdir = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/RotterdamStudy/MarjoleinVisit3/CellTypeSpecificityTestOutputMT2/";
//
//            getColumn(query, outdir, matrixIn);

            String matrixIn = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/2013-06-26-MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
            String query = "CellTypeInteractionZScore";
            String outdir = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/2013-06-26-MetaAnalysis/";

            getRow(query, outdir, matrixIn);
            
            
            matrixIn = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/2013-06-26-MetaAnalysis/MetaAnalysisZScoreMatrixSampleSize.txt";
            query = "CellTypeInteractionZScore";
            outdir = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/2013-06-26-MetaAnalysis/SampleSize-ˆ";

            getRow(query, outdir, matrixIn);


        } catch (IOException e) {

            e.printStackTrace();
        }
    }

    private static void getColumn(String query, String outdir, String matrixIn) throws IOException {
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(matrixIn);
        Integer id = ds.hashCols.get(query);
        if (id == null) {
            System.err.println("Query not in dataset..");
            System.exit(-1);
        }
        TextFile out = new TextFile(outdir + "Vector-" + query + "-2.txt", TextFile.W);
        out.writeln("\t" + query);
        for (int r = 0; r < ds.nrRows; r++) {
            out.writeln(ds.rowObjects.get(r) + "\t" + ds.rawData[r][id]);
        }
        out.close();
    }

    private static void getRow(String query, String outdir, String matrixIn) throws IOException {
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(matrixIn);
        Integer id = ds.hashRows.get(query);
        if (id == null) {
            System.err.println("Query not in dataset..");
            System.exit(-1);
        }

        TextFile out = new TextFile(outdir + "Vector-" + query + "-2.txt", TextFile.W);
        out.writeln("\t" + query);
        for (int r = 0; r < ds.nrCols; r++) {
            out.writeln(ds.colObjects.get(r) + "\t" + ds.rawData[id][r]);
        }
        out.close();
    }
}
