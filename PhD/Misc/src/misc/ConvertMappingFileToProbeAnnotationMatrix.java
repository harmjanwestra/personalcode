/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ConvertMappingFileToProbeAnnotationMatrix {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            ConvertMappingFileToProbeAnnotationMatrix c = new ConvertMappingFileToProbeAnnotationMatrix();
            String fastafile = "/Volumes/iSnackHD/Data/ProbeAnnotation/IlluminaProbeAnnotation/AllProbesFasta/AllIlluminaProbes.fa.txt";
            String platformAnnotationFileDir = "/Volumes/iSnackHD/Data/ProbeAnnotation/IlluminaProbeAnnotation/HJ/";
            String[] platformAnnotationFile = new String[]{
                "HT12v3.txt",
                "HT12v4.txt",
                "HT12v4WGDASL.txt",
                "HumanRef-8.txt",
                "HumanRef-8v2.txt",
                "HumanRef-8v3.txt",
                "HumanWG-6.txt",
                "HumanWG-6v2.txt",
                "HumanWG-6v3.txt"
            };

            String probeMappingFile = "/Volumes/iSnackHD/Data/ProbeAnnotation/IlluminaProbeAnnotation/Annotation/all.intersect.paf_unique_combined.txt";
            String output = "/Volumes/iSnackHD/Data/ProbeAnnotation/IlluminaProbeAnnotation/2013-07-08-ProbeAnnotationSummary.txt";
            c.concat(fastafile, platformAnnotationFileDir, platformAnnotationFile, probeMappingFile, output);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void concat(String fastaFile, String platformAnnotationFileDir, String[] platformAnnotationFile, String probeMappingFile, String output) throws IOException {
        HashMap<String, String> sequenceToProbe = new HashMap<String, String>();
        TextFile tf = new TextFile(fastaFile, TextFile.R);
        String ln = tf.readLine();
        ArrayList<String> sequences = new ArrayList<String>();
        while (ln != null) {
            if (ln.startsWith(">")) {
                // annotation in this line
                String probe = ln.replaceAll(">", "");
                String seq = tf.readLine().toUpperCase();
                sequences.add(seq);
                sequenceToProbe.put(seq, probe);
            }
            ln = tf.readLine();
        }
        tf.close();


        ArrayList<HashMap<String, ArrayList<String>>> sequencesPerPlatform = new ArrayList<HashMap<String, ArrayList<String>>>();
        for (int i = 0; i < platformAnnotationFile.length; i++) {

            HashMap<String, ArrayList<String>> data = new HashMap<String, ArrayList<String>>();
            TextFile tf2 = new TextFile(platformAnnotationFileDir + platformAnnotationFile[i], TextFile.R);
            String[] elems = tf2.readLineElems(TextFile.tab);
            while (elems != null) {
                String seq = elems[1].toUpperCase();
                String probe = elems[0];
                ArrayList<String> probes = data.get(seq);
                if (probes == null) {
                    probes = new ArrayList<String>();

                }
                probes.add(probe);
                data.put(seq, probes);
                elems = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();
            sequencesPerPlatform.add(data);
        }


        HashMap<String, String> probeToAnnotation = new HashMap<String, String>();
        TextFile annotationFile = new TextFile(probeMappingFile, TextFile.R);
        String[] elems = annotationFile.readLineElems(TextFile.tab);
        while (elems != null) {
            String probe = elems[1];
            String annotation = elems[3] + "\t" + elems[4] + "-"+elems[5]+"\t-";
            probeToAnnotation.put(probe, annotation);
            elems = annotationFile.readLineElems(TextFile.tab);
        }
        annotationFile.close();

        TextFile tfOut = new TextFile(output, TextFile.W);
        String header = "Probe\tSequence\tChromosome\tChromosome Position\tGene";

        for (int i = 0; i < platformAnnotationFile.length; i++) {
            header += "\t" + platformAnnotationFile[i];
        }
        tfOut.writeln(header);
        for (String sequence : sequences) {
            String probe = sequenceToProbe.get(sequence);
            String annotation = probeToAnnotation.get(probe);
            if (annotation == null) {
                annotation = "-\t-\t-";
            }
            String metaProbeOutput = probe + "\t" + sequence + "\t" + annotation;
            for (int i = 0; i < platformAnnotationFile.length; i++) {
                HashMap<String, ArrayList<String>> probemap = sequencesPerPlatform.get(i);
                ArrayList<String> probes = probemap.get(sequence);
                if (probes == null) {
                    metaProbeOutput += "\t-";
                } else {
                    metaProbeOutput += "\t" + Strings.concat(probes.toArray(new String[0]), Strings.semicolon);
                }
            }
            tfOut.writeln(metaProbeOutput);
        }
        tfOut.close();




    }
}
