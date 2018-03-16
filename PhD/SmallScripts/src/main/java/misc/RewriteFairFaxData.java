/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class RewriteFairFaxData {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String in = "/Volumes/iSnackHD/AeroFS/FairFaxData/Non_background_corrected_Bcell_monocyte_arrays.txt.gz";
        String out = "/Volumes/iSnackHD/AeroFS/FairFaxData/Non_background_corrected_Bcell_monocyte_arrays-AVGSignalOnly.gz";
        
        String annotationfile = "/Volumes/iSnackHD/AeroFS/FairFaxData/bcell_monocyte_covariates.txt";
        String monocyteout = "/Volumes/iSnackHD/AeroFS/FairFaxData/Non_background_corrected_monocyte_arrays-AVGSignalOnly.gz";
        String bcellcyteout = "/Volumes/iSnackHD/AeroFS/FairFaxData/Non_background_corrected_Bcell_arrays-AVGSignalOnly.gz";
        
        try {
            RewriteFairFaxData f = new RewriteFairFaxData();
            //f.run(in, out);
//            f.countelems(out);
            f.split(out, annotationfile, monocyteout, bcellcyteout);

//            in = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/BenFairfaxGenotypeData/Monocytes/Monocyte_raw_relabelled.txt";
//            out = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/BenFairfaxGenotypeData/Monocytes/Monocyte_raw_relabelled-AVGSignalOnly.txt";
//            f.run(in, out);
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    public void run(String infile, String outfile) throws IOException {
        TextFile in = new TextFile(infile, TextFile.R);
        TextFile out = new TextFile(outfile, TextFile.W);
        Pattern space = Pattern.compile("\\s");

        System.out.println(in.readLine());
        System.out.println(in.readLine());
        System.out.println(in.readLine());
        System.out.println(in.readLine());
        System.out.println(in.readLine());
        System.out.println(in.readLine());
        System.out.println(in.readLine());
        
        String header = in.readLine();
        String[] headerElems = space.split(header);
        out.append(headerElems[1]);
        for (int i = 2; i < headerElems.length; i++) {
            String headElem = headerElems[i];
            if (headElem.contains("AVG_Signal")) {
                headElem = headElem.replaceAll("AVG_Signal-", "");
                headElem = headElem.replaceAll("AVG_Signal", "");
                headElem = headElem.replaceAll("Bcell", "");
                headElem = headElem.replaceAll(".CD_", "");

                out.append("\t");
                out.append(headElem);
            }
        }
        out.append("\n");

        String ln = in.readLine();
        while (ln != null) {
            String[] elems = space.split(ln);
            out.append(elems[1]);
            for (int i = 2; i < elems.length; i++) {
                if (headerElems[i].contains("AVG")) {
                    out.append("\t");
                    out.append(elems[i]);
                }
            }
            out.append("\n");
            ln = in.readLine();
        }
        in.close();
        out.close();
    }

    private void countelems(String out) throws IOException {
        TextFile tf = new TextFile(out, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        System.out.println(elems.length);
        tf.close();
    }

    private void split(String out, String annotationfile, String monocyteout, String bcellcyteout) throws IOException {
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(out);
        HashMap<String, String> arrayToClass = new HashMap<String, String>();
        HashMap<String, String> arrayToId = new HashMap<String, String>();
        
        TextFile tf = new TextFile(annotationfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        int nrMonocyte = 0;
        int nrBCell = 0;
        while(elems!=null){
            String array = elems[1];
            String id = elems[2];
            String sampletype = elems[3];
            
            if(sampletype.equals("monocyte")){
                nrMonocyte++;
            } else if(sampletype.equals("Bcells")){
                nrBCell++;
            }
            
            arrayToClass.put(array, sampletype);
            arrayToId.put(array, id);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        
        System.out.println("Nr B:"+nrBCell);
        System.out.println("Nr M:"+nrMonocyte);
        
        ds.transposeDataset();
        double[][] datamonocyte = new double[nrMonocyte][ds.nrCols];
        ArrayList<String> monocyteSample = new ArrayList<String>();
        double[][] databcell = new double[nrBCell][ds.nrCols];
        ArrayList<String> bcellSample = new ArrayList<String>();
        for(int r=0;r<ds.nrRows;r++){
            String sample = ds.rowObjects.get(r);
            String id = arrayToId.get(sample);
            String sampletype = arrayToClass.get(sample);
            if(sampletype.equals("monocyte")){
                datamonocyte[monocyteSample.size()] = ds.rawData[r];
                monocyteSample.add(id);
            } else if(sampletype.equals("Bcells")){
                databcell[bcellSample.size()] = ds.rawData[r];
                bcellSample.add(id);
            }
        }
        
        DoubleMatrixDataset<String, String> outfile = new DoubleMatrixDataset<String, String>();
        outfile.rawData = datamonocyte;
        outfile.rowObjects = monocyteSample;
        outfile.colObjects = ds.colObjects;
        outfile.recalculateHashMaps();
        outfile.transposeDataset();
        outfile.save(monocyteout);
        
        outfile.rawData = databcell;
        outfile.rowObjects = bcellSample;
        outfile.colObjects = ds.colObjects;
        outfile.recalculateHashMaps();
        outfile.transposeDataset();
        outfile.save(bcellcyteout);
        
        
    }
}
