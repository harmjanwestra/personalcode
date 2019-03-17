package nl.harmjanwestra.playground;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class VCFRSFilter {


    public static void main(String[] args) {


        if (args.length < 3) {
            System.out.println("Usage: invcf rsids outvcf");
        } else {


            try {
                System.out.println("Reading: " + args[1]);
                TextFile tf = new TextFile(args[1], TextFile.R);
                HashSet<String> query = new HashSet<String>();
                query.addAll(tf.readAsArrayList());
                tf.close();
                System.out.println(query.size() + " SNPs to be selected");


                System.out.println("Input: " + args[0]);
                System.out.println("Output: " + args[2]);
                TextFile tf2 = new TextFile(args[0], TextFile.R);
                TextFile tf3 = new TextFile(args[2], TextFile.W);
                String ln = tf2.readLine();
                int lnctr = 0;
                int written = 0;
                while (ln != null) {

                    if (ln.startsWith("#")) {
                        tf3.writeln(ln);
                        written++;
                    } else {
                        String[] elems = Strings.subsplit(ln, Strings.tab, 0, 3);
                        String rs = elems[2];
                        if (query.contains(rs)) {
                            tf3.writeln(ln);
                            written++;
                        }
                    }
                    ln = tf2.readLine();
                    lnctr++;
                    if (lnctr % 1000000 == 0) {
                        System.out.print("\r" + lnctr + " lines parsed sofar.");
                    }
                }
                System.out.println();
                System.out.println("Done. " + lnctr + " lines parsed in total, " + written + " lines written. ");
                tf2.close();
                tf3.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
