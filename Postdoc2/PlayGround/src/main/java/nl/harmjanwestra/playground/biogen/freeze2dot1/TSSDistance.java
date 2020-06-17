package nl.harmjanwestra.playground.biogen.freeze2dot1;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

public class TSSDistance {

    public static void main(String[] args) {
        String efilefolder = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cis\\";

        String[] names = new String[]{
                "2020-05-26-Basalganglia-EUR-IterationITER-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "2020-05-26-Cerebellum-EUR-IterationITER-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "2020-05-26-Cortex-AFR-IterationITER-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "2020-05-26-Cortex-EAS-IterationITER-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "2020-05-26-Cortex-EUR-IterationITER-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "2020-05-26-Hippocampus-EUR-IterationITER-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "2020-05-26-Spinalcord-EUR-IterationITER-eQTLProbesFDR0.05-ProbeLevel.txt.gz"
        };

        TSSDistance d = new TSSDistance();
        try {
            int range = 50000;
            String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\tssdistance\\test" + range + ".txt";
            d.run(efilefolder, names, 4, range, out);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String efiledir, String[] names, int maxiter, int range, String out) throws IOException {

        ArrayList<String> files = new ArrayList<String>();
        ArrayList<String> dsnames = new ArrayList<String>();
        for (int i = 0; i < names.length; i++) {
            for (int j = 1; j < maxiter; j++) {
                String file = efiledir + names[i].replaceAll("ITER", "" + j);
                System.out.println("Checking: " + file);
                if (Gpio.exists(file)) {
                    files.add(file);
                    dsnames.add(names[i].replaceAll("ITER", "" + j));
                }
            }
        }


        System.out.println(files.size() + " files");
        System.out.println();
        int nrbins = 40;
        double[][] bins = new double[nrbins + 1][files.size()];

        for (int f = 0; f < files.size(); f++) {
            TextFile tf = new TextFile(files.get(f), TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String snppos = elems[3];
                String genepos = elems[6];
                Integer snp = Integer.parseInt(snppos);
                Integer gene = Integer.parseInt(genepos);

                int delta = snp - gene;
                double bind = ((double) (delta + (range / 2)) / range) * nrbins;
                if (bind < 0) {
                    bind = 0;
                }
                if (bind >= nrbins) {
                    bind = nrbins;
                }
                int bini = (int) Math.floor(bind);
//                System.out.println(snp + "\t" + gene + "\t" + delta + "\t" + bini);
                bins[bini][f]++;
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
        }

        // convert to percentages
        for (int f = 0; f < files.size(); f++) {
            double sum = 0;
            for (int i = 0; i < bins.length; i++) {
                sum += bins[i][f];
            }
            for (int i = 0; i < bins.length; i++) {
                bins[i][f] /= sum;
                if (Double.isNaN(bins[i][f])) {
                    bins[i][f] = 0;
                }
            }
        }


        TextFile tfo = new TextFile(out, TextFile.W);
        String header = "Bin\tbp\t" + Strings.concat(dsnames, Strings.tab);
        tfo.writeln(header);
        int bpPerBin = range / nrbins;
        for (int i = 0; i < bins.length; i++) {
            int bp = (bpPerBin * i) - (range / 2);
            tfo.writeln(i + "\t" + bp + "\t" + Strings.concat(bins[i], Strings.tab));
            System.out.println(i + "\t" + bp + "\t" + Strings.concat(bins[i], Strings.tab));
        }
        tfo.close();
    }
}
