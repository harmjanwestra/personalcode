/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import eqtlmappingpipeline.normalization.Normalizer;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ConvertArrayExpressDataToExpressionData {

    HashMap<String, String> ilmnToArrayAddress = new HashMap<String, String>();

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        try {
            ConvertArrayExpressDataToExpressionData d = new ConvertArrayExpressDataToExpressionData();
            d.loadProbeAnnotation("/Volumes/iSnackHD/Data/ProbeAnnotation/ProbeAnnotationAccordingToIllumina/HumanWG-6_V2_0_R4_11223189_A.txt", "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/WG6ProbeAnnotation.txt");

//            d.run1("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpressionDataRaw/CHB.raw_data.txt.gz",
//                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataCHB/ExpressionDataWithDuplicates.txt.gz");
//
//            d.run1("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpressionDataRaw/GIH.raw_data.txt.gz",
//                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataGIH/ExpressionDataWithDuplicates.txt.gz");
//
//            d.run1("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpressionDataRaw/JPT.raw_data.txt.gz",
//                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataJPT/ExpressionDataWithDuplicates.txt.gz");
//
//            d.run1("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpressionDataRaw/LWK.raw_data.txt.gz",
//                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataLWK/ExpressionDataWithDuplicates.txt.gz");
//
//            d.run1("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpressionDataRaw/MEX.raw_data.txt.gz",
//                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataMEX/ExpressionDataWithDuplicates.txt.gz");
//
//            d.run1("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpressionDataRaw/MKK.raw_data.txt.gz",
//                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataMKK/ExpressionDataWithDuplicates.txt.gz");
//
//            d.run1("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpressionDataRaw/YRI.raw_data.txt.gz",
//                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataYRI/ExpressionDataWithDuplicates.txt.gz");

            // combine
            
            d.run("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap-Stranger/RAW/ExpressionData.txt",
                  "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap-Stranger/RAW/ExpressionDataWithoutDuplicates.txt.gz");
//            
//            d.run("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataCHB/ExpressionDataWithDuplicates.txt.gz",
//                  "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataCHB/ExpressionDataWithOutDuplicatesWG6.txt.gz");
//
//            d.run("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataGIH/ExpressionDataWithDuplicates.txt.gz",
//                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataGIH/ExpressionDataWithOutDuplicatesWG6.txt.gz");
//
//            d.run("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataJPT/ExpressionDataWithDuplicates.txt.gz",
//                  "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataJPT/ExpressionDataWithOutDuplicatesWG6.txt.gz");
//
//            d.run("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataLWK/ExpressionDataWithDuplicates.txt.gz",
//                  "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataLWK/ExpressionDataWithOutDuplicatesWG6.txt.gz");
//
//            d.run("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataMEX/ExpressionDataWithDuplicates.txt.gz",
//                  "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataMEX/ExpressionDataWithOutDuplicatesWG6.txt.gz");
//
//            d.run("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataMKK/ExpressionDataWithDuplicates.txt.gz",
//                  "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataMKK/ExpressionDataWithOutDuplicatesWG6.txt.gz");
//
//            d.run("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataYRI/ExpressionDataWithDuplicates.txt.gz",
//                  "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataYRI/ExpressionDataWithOutDuplicatesWG6.txt.gz");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run1(String in, String out) throws IOException {
        
        
        
        TextFile inFile = new TextFile(in, TextFile.R);
        String[] headerElems = inFile.readLineElems(TextFile.tab);

        boolean[] colToInclude = new boolean[headerElems.length];
        ArrayList<String> sampleNames = new ArrayList<String>();
        HashSet<String> uniqueNames = new HashSet<String>();
        HashSet<String> uniqueSamples = new HashSet<String>();
        HashMap<String, Integer> firstSampleId = new HashMap<String, Integer>();

        for (int i = 1; i < colToInclude.length; i++) {
            if (headerElems[i].endsWith(".AVG_Signal")) {
                colToInclude[i] = true;
                String name = headerElems[i];

                if (!uniqueNames.contains(name)) {
                    sampleNames.add(name);
                    firstSampleId.put(name, uniqueNames.size());
                    uniqueNames.add(name);
                    uniqueSamples.add(name);

                } else {
//                    System.out.print("Duplicate sample: "+name+". ");
//                    while (uniqueNames.contains(name)) {
//                        name+="I";
//                    }
//                    uniqueNames.add(name);
//                    System.out.print("Replaced with: "+name+"\n");
//                    sampleNames.add(name);
                }
            }
        }

        System.out.println("Total samples: " + sampleNames.size() + "\tUnique samples: " + uniqueSamples.size());

        TextFile outFile = new TextFile(out, TextFile.W);

        String headerOut = Strings.concat(sampleNames, Strings.tab);
        outFile.writeln("Probe\t" + headerOut);

        String[] elems = inFile.readLineElems(TextFile.tab);
        while (elems != null) {

            double[] vals = new double[uniqueSamples.size()];


            for (int i = 1; i < elems.length; i++) {
                if (colToInclude[i]) {
                    String name = headerElems[i];
                    Integer id = firstSampleId.get(name);
                    vals[id] += Double.parseDouble(elems[i]);
                }
            }

            String[] valsOut = new String[vals.length];
            for (int i = 0; i < vals.length; i++) {
//                vals[i] /= 2;
                valsOut[i] = "" + vals[i];
            }

            String outStr = Strings.concat(valsOut, Strings.tab);
            outFile.writeln(elems[0] + "\t" + outStr);

            elems = inFile.readLineElems(TextFile.tab);
        }

        outFile.close();
        inFile.close();


    }

    public void run(String in, String out) throws IOException {
        
        Normalizer m = new Normalizer();
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String,String>(in);
        in = m.quantileNormalize(ds, in);
        in+=".txt.gz";
        TextFile inFile = new TextFile(in, TextFile.R);
        String[] headerElems = inFile.readLineElems(TextFile.tab);

        boolean[] colToInclude = new boolean[headerElems.length];
        ArrayList<String> sampleNames = new ArrayList<String>();
        HashSet<String> uniqueNames = new HashSet<String>();
        HashSet<String> uniqueSamples = new HashSet<String>();
        HashMap<String, Integer> firstSampleId = new HashMap<String, Integer>();

        for (int i = 1; i < colToInclude.length; i++) {
           // if (headerElems[i].endsWith(".AVG_Signal")) {
                colToInclude[i] = true;
                String name = headerElems[i].split("_")[0];

                if (!uniqueNames.contains(name)) {
                    sampleNames.add(name);
                    firstSampleId.put(name, uniqueNames.size());
                    uniqueNames.add(name);
                    uniqueSamples.add(name);

                } else {
//                    System.out.print("Duplicate sample: "+name+". ");
//                    while (uniqueNames.contains(name)) {
//                        name+="I";
//                    }
//                    uniqueNames.add(name);
//                    System.out.print("Replaced with: "+name+"\n");
//                    sampleNames.add(name);
                }
            //}
        }

        System.out.println("Total samples: " + sampleNames.size() + "\tUnique samples: " + uniqueSamples.size());

        TextFile outFile = new TextFile(out, TextFile.W);

        String headerOut = Strings.concat(sampleNames, Strings.tab);
        outFile.writeln("Probe\t" + headerOut);

        String[] elems = inFile.readLineElems(TextFile.tab);
        while (elems != null) {

            double[] vals = new double[uniqueSamples.size()];


            for (int i = 1; i < elems.length; i++) {
                if (colToInclude[i]) {
                    String name = headerElems[i].split("_")[0];
                    Integer id = firstSampleId.get(name);
                    vals[id] += Double.parseDouble(elems[i]);
                }
            }

            String[] valsOut = new String[vals.length];
            for (int i = 0; i < vals.length; i++) {
                vals[i] /= 2;
                valsOut[i] = "" + vals[i];
            }

            String outStr = Strings.concat(valsOut, Strings.tab);
            outFile.writeln(elems[0] + "\t" + outStr);

            elems = inFile.readLineElems(TextFile.tab);
        }

        outFile.close();
        inFile.close();


    }

    private void loadProbeAnnotation(String string, String string0) throws IOException {
        TextFile in = new TextFile(string, TextFile.R);
//        TextFile out = new TextFile(string0, TextFile.R);

        String[] elems = in.readLineElems(TextFile.tab);
        int pbCol = -1;
        int arrCol = -1;
        for (int i = 0; i < elems.length; i++) {
            if (elems[i].equals("Probe_Id")) {
                pbCol = i;
            }
            if (elems[i].equals("Array_Address_Id")) {
                arrCol = i;
            }
        }

        System.out.println("Probe: " + pbCol);
        System.out.println("ArrayID: " + arrCol);
        int nrProc = 1;
        elems = in.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length > pbCol && elems.length > arrCol) {

                String arr = elems[arrCol];
                while (arr.startsWith("0")) {
                    arr = arr.substring(1);
                }

                if (ilmnToArrayAddress.containsKey(elems[pbCol]) || ilmnToArrayAddress.containsValue(arr)) {
                    System.out.println("Duplicate probe!" + elems[pbCol] + "\t" + arr);
                } else {
                    ilmnToArrayAddress.put(elems[pbCol], arr);
                }

            }
            elems = in.readLineElems(TextFile.tab);
            if (nrProc % 1000 == 0) {
                System.out.println(nrProc);
            }
            nrProc++;
        }

        in.close();
    }
}
