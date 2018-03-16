/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.NavigableMap;
import java.util.TreeMap;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harm-jan
 */
public class FraminghamDataIntegrator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String[] eqtlfiles = new String[]{"D:\\AeroFS\\Framingham\\EQTLListBest100NonRedundant.txt", "D:\\AeroFS\\Framingham\\EQTLList.txt"};
        String probeAnnotationFile = "D:\\AeroFS\\cellTypeeQTL\\DataFiles\\2013-07-18-ProbeAnnotationFile.txt";
        String snpAnnotationFile = "D:\\AeroFS\\Framingham\\2013-08-05-SNPMappings_Filtered.txt.gz";
        String file = "/Volumes/iSnackHD/Data/Projects/Framingham/TransEQTLReplication/EQTLListBest100NonRedundant.txt-Appended.txt";
        String outfile = "/Volumes/iSnackHD/Data/Projects/Framingham/TransEQTLReplication/EQTLListBest100NonRedundant.txt-Appended-SNPProbeCombos.txt";

        try {
            FraminghamDataIntegrator f = new FraminghamDataIntegrator();
//            f.convert(eqtlfiles, probeAnnotationFile, snpAnnotationFile);
            f.split(file, outfile);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void convert(String[] input, String probeannotationfile, String snpannotationfile) throws IOException {
        // load file to be anotated..
        for (String file : input) {
            TextFile tffile = new TextFile(file, TextFile.R);
            ArrayList<String> alldata = new ArrayList<String>();
            String header = tffile.readLine();
            String line = tffile.readLine();
            while (line != null) {

                alldata.add(line);

                line = tffile.readLine();
            }
            tffile.close();

            String[] lines = alldata.toArray(new String[0]);

            for (byte chr = 0; chr < 26; chr++) {

                System.out.println(chr);
                NavigableMap<Integer, String> snps = new TreeMap<Integer, String>();

                // load SNP annotation
                TextFile tf = new TextFile(snpannotationfile, TextFile.R);
                String[] elems = tf.readLineElems(TextFile.tab);
                int ctr = 0;
                while (elems != null) {

                    byte chr2 = ChrAnnotation.parseChr(elems[0]);
                    Integer chrPos = Integer.parseInt(elems[1]);
                    String rs = new String("" + elems[2]);
                    if (chr2 == chr) {
                        snps.put(chrPos, rs);
                    }
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();



                for (int i = 0; i < lines.length; i++) {
                    line = lines[i];
                    String[] lineElems = line.split("\t");
                    byte chr2 = ChrAnnotation.parseChr(lineElems[1]);
                    if (chr2 == chr) {
                        Integer chrPos = Integer.parseInt(lineElems[2]);

                        byte transcriptchr = ChrAnnotation.parseChr(lineElems[3]);
                        Integer transcriptStartPos = Integer.parseInt(lineElems[4]);
                        Integer transcriptStopPos = Integer.parseInt(lineElems[5]);
                        String rs = snps.get(chrPos);
                        if (rs != null) {
                            lines[i] += "\t" + rs;
                        } else {
                            lines[i] += "\t" + null;
                        }

                        TextFile probefile = new TextFile(probeannotationfile, TextFile.R);
                        probefile.readLine();
                        String[] probeElems = probefile.readLineElems(TextFile.tab);
                        ArrayList<String> probesList = new ArrayList<String>();
                        while (probeElems != null) {
                            byte probechr = ChrAnnotation.parseChr(probeElems[2]);
                            if (transcriptchr == probechr) {
                                String pos = probeElems[3];
                                String[] posElems = pos.split("-");
                                String start = posElems[0].split(":")[0];
                                Integer startI = Integer.parseInt(start);
                                String[] stopelems = posElems[posElems.length - 1].split(":");
                                Integer stopI = Integer.parseInt(stopelems[stopelems.length - 1]);
                                if (startI >= transcriptStartPos && stopI <= transcriptStopPos) {
                                    probesList.add(probeElems[5]);
                                }
                            }
                            probeElems = probefile.readLineElems(TextFile.tab);
                        }
                        probefile.close();

                        lines[i] += "\t" + Strings.concat(probesList, Strings.comma);
                    }
                }
            }

            TextFile output = new TextFile(file + "-Appended.txt", TextFile.W);
            output.writeln(header);
            for (String s : lines) {
                output.writeln(s);
            }
            output.close();


        }


    }

    private void split(String file, String outfile) throws IOException {
        TextFile tf = new TextFile(file, TextFile.R);
        TextFile tfout = new TextFile(outfile, TextFile.W);
        String[] elems = tf.readLineElems(TextFile.tab);
        HashSet<Pair<String, String>> pairs = new HashSet<Pair<String, String>>();
        while (elems != null) {
            if (elems.length < 8) {
                System.out.println(Strings.concat(elems, Strings.tab));
            } else {
                String snp = elems[6];
                String probes = elems[7];
                String[] probeElems = probes.split(",");

                for (String probe : probeElems) {

                    if (!probe.equals("-")) {
                        Pair<String, String> p = new Pair<String, String>(snp, probe);
                        if (!pairs.contains(p)) {
                            tfout.writeln(snp + "\t" + probe);
                            pairs.add(p);
                        }
                    }
                }
            }

            elems = tf.readLineElems(TextFile.tab);
        }
        tfout.close();
        tf.close();
    }
}
