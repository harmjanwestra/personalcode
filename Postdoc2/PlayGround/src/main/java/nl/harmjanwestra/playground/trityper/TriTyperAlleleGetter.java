package nl.harmjanwestra.playground.trityper;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.IOException;

public class TriTyperAlleleGetter {

    public static void main(String[] args) {
        TriTyperAlleleGetter g = new TriTyperAlleleGetter();
        if (args.length < 2) {
            System.out.println("Usage: trityperdir outfile");
        } else {
            try {
                g.run(args[0], args[1]);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public void run(String indir, String out) throws IOException {
        System.out.println("TriTyper to bed converter");
        System.out.println("in " + indir);
        System.out.println("out: " + out);
        TriTyperGenotypeData ds = new TriTyperGenotypeData(indir);
        String[] snps = ds.getSNPs();
        SNPLoader loader = ds.createSNPLoader(1024);



        ProgressBar pb = new ProgressBar(snps.length, "Convering alleles..");
        TextFile tf = new TextFile(out, TextFile.W);
        for (int s = 0; s < snps.length; s++) {
            SNP snp = ds.getSNPObject(s);
            loader.loadGenotypes(snp);
            byte[] allelesb = snp.getAlleles();
            String[] alleles = new String[allelesb.length];
            for (int b = 0; b < alleles.length; b++) {
                alleles[b] = BaseAnnot.toString(allelesb[b]);
            }
            tf.writeln("chr" + snp.getChr() + "\t" + snp.getChrPos() + "\t" + (snp.getChrPos() + 1) + "\t" + snp.getName() + "\t" + Strings.concat(alleles, Strings.tab));

            snp.clearGenotypes();
            pb.iterate();
        }
        pb.close();
        tf.close();

    }

}
