package nl.harmjanwestra.playground.biogen.gwasfilter;

import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class CombineAndConquer {

    public static void main(String[] args) {
        String gwas = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gwas_catalog_v1.0.2-associations_e98_r2020-02-25.tsv.gz";
        String efile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\trans-cortex-EURandAFR\\eQTLsFDR0.05.txt.gz";
        String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\trans-cortex-EURandAFR\\eQTLsFDR0.05-traits.txt.gz";
        CombineAndConquer c = new CombineAndConquer();
        try {
            c.run(gwas, efile, out);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    public void run(String gwas, String efile, String out) throws IOException {
        GWASCatalog ctl = new GWASCatalog(gwas);


        HashMap<String, ArrayList<EQTL>> snpToEQTL = new HashMap<>();
        QTLTextFile tf = new QTLTextFile(efile, TextFile.R);
        EQTL[] set = tf.read();
        tf.close();

        for (EQTL e : set) {
            String rsName = e.getRsName().split(":")[2];
            ArrayList<EQTL> list = snpToEQTL.get(rsName);
            if (list == null) {
                list = new ArrayList<>();
            }

            list.add(e);
            snpToEQTL.put(rsName, list);

        }

        TextFile outf = new TextFile(out, TextFile.W);

        GWASSNP[] allsnps = ctl.getSnpsArray();
        for (GWASSNP snp : allsnps) {
            String name = snp.getName();
            ArrayList<EQTL> elist = snpToEQTL.get(name);
            if (elist != null) {
                GWASTrait[] diseases = snp.getAssociatedTraitsArray();
                for (EQTL e : elist) {
                    for (GWASTrait t : diseases) {

                        String ln = name + "\t" + e.getProbe() + "\t" + e.getProbeHUGO() + "\t" + e.getFDR() + "\t" + t.getCleanName() + "\t" + t.getName();
                        outf.writeln(ln);

                    }
                }
            }

        }

        outf.close();

    }
}
