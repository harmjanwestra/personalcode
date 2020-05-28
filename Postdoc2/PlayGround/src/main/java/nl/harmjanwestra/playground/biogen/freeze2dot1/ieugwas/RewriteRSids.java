package nl.harmjanwestra.playground.biogen.freeze2dot1.ieugwas;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class RewriteRSids {

    public static void main(String[] args) {
        String allassocfile = "U:\\IEUGWAS\\2020-05-03-allTopAssociations-wgwascatalog-wALS.txt.gz";
        String allassocfileo = "U:\\IEUGWAS\\2020-05-03-allTopAssociations-wgwascatalog-wALS-MetaBrain2dot1IDs.txt.gz";
        String topassocfile = "U:\\IEUGWAS\\2020-05-03-topHitsUniqueRSids-wgwascatalog-wALS.txt.gz";
        String topassocfileo = "U:\\IEUGWAS\\2020-05-03-topHitsUniqueRSids-wgwascatalog-wALS-MetaBrain2dot1IDs.txt.gz";
        String refrs = "U:\\IEUGWAS\\allSNPs.txt.gz";
        RewriteRSids r = new RewriteRSids();
        try {
            r.rewrite(refrs, allassocfile, allassocfileo, 1);
            r.rewrite(refrs, topassocfile, topassocfileo, 0);
        } catch (IOException e) {
            e.printStackTrace();
        }



    }

    public void rewrite(String refrs, String ieugwasfile, String ieugwasfileout, int col) throws IOException {

        HashMap<String, String> rstors = new HashMap<String, String>();
        TextFile tf = new TextFile(refrs, TextFile.R);
        String ln = tf.readLine();
        while (ln != null) {
            String[] elems = ln.split(":");
            String rs = elems[2];
            rstors.put(rs, ln);
            ln = tf.readLine();
        }
        tf.close();

        TextFile in = new TextFile(ieugwasfile, TextFile.R);
        TextFile outf = new TextFile(ieugwasfileout, TextFile.W);
        outf.writeln(in.readLine());
        String[] elems = in.readLineElems(TextFile.tab);
        while (elems != null) {
            String id = elems[col];
            String rs = rstors.get(id);
            if (rs != null) {
                elems[col] = rs;
                outf.writeln(Strings.concat(elems, Strings.tab));
            }
            elems = in.readLineElems(TextFile.tab);
        }
        outf.close();
        in.close();

    }


}
