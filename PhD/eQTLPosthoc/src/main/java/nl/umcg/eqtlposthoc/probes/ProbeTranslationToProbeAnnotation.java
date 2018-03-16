/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.probes;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ProbeTranslationToProbeAnnotation {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String probeTranslationFile = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-ProbeAnnotationFile.txt";
        String probeAnnotationDir = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-";
        ProbeTranslationToProbeAnnotation pb = new ProbeTranslationToProbeAnnotation();



        try {
            pb.run(probeTranslationFile, probeAnnotationDir);

            String oldFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-H8v2.txt";
            String newFile = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-HumanRef-8v2.txt";
            String output = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-HumanRef-8v2-HT12v3Rewrite.txt";

            pb.rewriteH8v2Annotation(oldFile, newFile, output);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String probeTranslationFile, String probeAnnotationDir) throws IOException {
        TextFile tf = new TextFile(probeTranslationFile, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        for (int i = 5; i < header.length; i++) {
            TextFile tfout = new TextFile(probeAnnotationDir + header[i], TextFile.W);
            tf.close();
            tf.open();
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            String platform = header[i].replace(".txt", "");
            while (elems != null) {
                String probeStr = elems[i];
                if (!probeStr.equals("-")) {

                    String[] probeElems = probeStr.split(";");
                    for (String probe : probeElems) {

                        String metaprobe = elems[0];
                        String seq = elems[1];
                        String chr = elems[2];
                        String location = elems[3];
                        String gene = elems[4];

                        String chrStart = "-";
                        String chrEnd = "-";
                        if (!location.equals("-1") && !location.equals("-")) {
                            String[] posElems = location.split(":");
                            String[] posElems2 = posElems[0].split("-");
                            String[] posElems3 = posElems[posElems.length - 1].split("-");

                            chrStart = posElems2[0];
                            chrEnd = posElems3[posElems3.length - 1];
                        }

                        //Platform	HT12v3-ArrayAddress	Symbol	Chr	ChrStart	ChrEnd	Probe	Seq
                        tfout.writeln(platform + "\t" + probe + "\t" + gene + "\t" + chr + "\t" + chrStart + "\t" + chrEnd + "\t" + metaprobe + "\t" + seq);
                    }
                }

                elems = tf.readLineElems(TextFile.tab);
            }
            tfout.close();
        }
        tf.close();
    }

    private void rewriteH8v2Annotation(String oldFile, String newFile, String output) throws IOException {
        HashMap<String, String> sequenceToProbeName = new HashMap<String, String>();

        TextFile tf = new TextFile(oldFile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String probe = elems[1];
            String seq = elems[elems.length - 1];
            sequenceToProbeName.put(seq, probe);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        TextFile out = new TextFile(output, TextFile.W);
        TextFile in = new TextFile(newFile, TextFile.R);

        elems = in.readLineElems(TextFile.tab);
        out.writeln(Strings.concat(elems, Strings.tab));
        while (elems != null) {
            String seq = elems[elems.length - 1];
            elems[1] = sequenceToProbeName.get(seq);
            if (elems[1] != null) {
                out.writeln(Strings.concat(elems, Strings.tab));
            }
            elems = in.readLineElems(TextFile.tab);
        }
        in.close();
        out.close();

    }
}
