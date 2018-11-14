package nl.harmjanwestra.playground.biogen.cleangeuvadis;

import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class CleanMat {

    public static void main(String[] args) {

        String europeans = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p35va-europeans.txt";
        String in = "D:\\geuvadis\\GD660.GeneQuantCount.txt.gz";
        String out = "D:\\geuvadis\\GD660.GeneQuantCount-EUR.txt.gz";

        CleanMat c = new CleanMat();
        try {
            c.run(in, europeans, out);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String in, String samplefilter, String out) throws IOException {

        HashSet<String> incSamples = null;
        if (samplefilter != null) {
            TextFile tf = new TextFile(samplefilter, TextFile.R);
            ArrayList<String> list = tf.readAsArrayList();
            tf.close();
            incSamples = new HashSet<>();
            incSamples.addAll(list);
        }

        System.out.println("parsing: " + in);

        TextFile tf3 = new TextFile(out + "GeneAnnotation.txt", TextFile.W);
        tf3.writeln("Platform\tArrayAddress\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tSeq");
        TextFile tf = new TextFile(in, TextFile.R);
        TextFile tf2 = new TextFile(out, TextFile.W);
        String header = tf.readLine();
        String[] headerElems = header.split("\t");

        boolean[] includecol = null;
        int nrsamplesoverlap = 0;
        if (incSamples == null) {
            tf2.writeln(header);
        } else {
            String newheader = "-";
            includecol = new boolean[headerElems.length];
            includecol[0] = true;
            nrsamplesoverlap = headerElems.length - 1;
            for (int s = 1; s < headerElems.length; s++) {

                String sampleid = headerElems[s].split("\\.")[0];
                headerElems[s] = sampleid;

                if (incSamples.contains(sampleid)) {
                    includecol[s] = true;
                    nrsamplesoverlap++;
                }
            }
            newheader = Strings.concat(headerElems, includecol, Strings.tab);
            tf2.writeln(newheader);
        }

//        tf.readLine(); //  unmapped
//        tf.readLine(); //  multimapped
//        tf.readLine(); // nofeature
//        tf.readLine(); // ambiguous
        String[] elems = tf.readLineElems(Strings.tab); // ensg.1
        int written = 0;
        int read = 0;
        while (elems != null) {
//				elems[0] = elems[0]; // .split("\\.")[0];
            int nrzeros = 0;
            for (int q = 1; q < elems.length; q++) {

                if (includecol != null && includecol[q]) {
                    Double d = Double.parseDouble(elems[q]);
                    if (d == 0d) {
                        nrzeros++;
                    }
                } else if (includecol == null) {
                    Double d = Double.parseDouble(elems[q]);
                    if (d == 0d) {
                        nrzeros++;
                    }
                }
            }
            if (nrzeros != nrsamplesoverlap) {
                tf3.writeln("1kg\t" + elems[0] + "\t" + elems[0] + "\t" + elems[2] + "\t" + elems[3] + "\t" + elems[3] + "\t" + elems[0] + "\t-");
                if (includecol != null) {
                    tf2.writeln(Strings.concat(elems, includecol, Strings.tab));
                } else {
                    tf2.writeln(Strings.concat(elems, Strings.tab));
                }

            }


            elems = tf.readLineElems(Strings.tab);
        }
        tf2.close();
        tf.close();
        tf3.close();

    }
}
