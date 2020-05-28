package nl.harmjanwestra.playground.biogen.freeze2dot1;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;

public class ConvertEQTLMappingPipelineZScoreMatrix {

    public static void main(String[] args) {

        String rowfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\trans-cortex-EURandAFR-zscorematrix\\MetaAnalysis-RowNames.txt.gz";
        String colfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\trans-cortex-EURandAFR-zscorematrix\\MetaAnalysis-ColNames.txt.gz";
        String matrix = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\trans-cortex-EURandAFR-zscorematrix\\MetaAnalysis.dat";
        String outputZ = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\trans-cortex-EURandAFR-zscorematrix\\MetaAnalysis-Z.txt.gz";
        String outputN = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\trans-cortex-EURandAFR-zscorematrix\\MetaAnalysis-N.txt.gz";

        ConvertEQTLMappingPipelineZScoreMatrix c = new ConvertEQTLMappingPipelineZScoreMatrix();
        try {
            c.convert(rowfile, colfile, matrix, outputZ, outputN);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    class SNP {
        String id;
        String alleles;
        String assessed;
        int n;
    }

    public void convert(String rowfile, String colfile, String matrix, String outputZ, String outputN) throws IOException {

        ArrayList<SNP> snps = new ArrayList<SNP>();
        TextFile tf = new TextFile(rowfile, TextFile.R);
        tf.readLine();//skip header
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            // 1:960326:rs13303327:G_A A/G     G       G       531     -       -       -       19373   -
            SNP snp = new SNP();
            snp.id = elems[0];
            snp.alleles = elems[1];
            snp.assessed = elems[3];
            snp.n = Integer.parseInt(elems[4]);
            snps.add(snp);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        tf = new TextFile(colfile, TextFile.R);
        ArrayList<String> genes = tf.readAsArrayList();
        tf.close();

        ArrayList<String> tmpgenes = new ArrayList<>();
        for (int g = 0; g < genes.size(); g++) {
            String name = genes.get(g).split("\\.")[0];
            tmpgenes.add(name);
        }
        genes = tmpgenes;

        System.out.println(genes.size());
        System.out.println(snps.size());
        // output format:
        // snp alleles alleleassessed gene1_gene1 gene2_gene2 etc

        // single int as header
        long expectedsize = 4 + (long) snps.size() * (long) genes.size() * (long) 4;
        File f = new File(matrix);
        long actuallength = f.length();
        if (expectedsize != actuallength) {
            System.err.println("error: was expecting: " + expectedsize + " but found " + actuallength + " diff: " + (expectedsize - actuallength) + " bytes difference");
            System.exit(-1);
        }
        BinaryFile bf = new BinaryFile(matrix, BinaryFile.R); // does this have a magic number or header?

        int intval = bf.readInt();
        System.out.println("Type of file: " + intval);

        TextFile tfout = new TextFile(outputZ, TextFile.W);
        TextFile tfout2 = new TextFile(outputN, TextFile.W);
        tfout.writeln("snp\talleles\talleleassessed\t" + Strings.concat(genes, Strings.tab));
        ProgressBar pb = new ProgressBar(snps.size());
        byte[] readbuffer = new byte[genes.size() * 4 * 100];
        byte[] snpbuffer = new byte[genes.size() * 4];
        for (int s = 0; s < snps.size(); s++) {

            if (s % 100 == 0) {
                if (s + 100 > snps.size()) {
                    int n = snps.size() - s;
                    readbuffer = new byte[genes.size() * 4 * n];
                }
                bf.read(readbuffer);
            }

            SNP snp = snps.get(s);
            String ln = snp.id + "\t" + snp.alleles + "\t" + snp.assessed;
            String ln2 = ln;

//System.arraycopy(bZs, offset, bytesToRead, 0, readlen);
            int offset = (s % 100) * genes.size() * 4;
            System.arraycopy(readbuffer, offset, snpbuffer, 0, snpbuffer.length);

            ByteBuffer bytebuffer = ByteBuffer.wrap(snpbuffer);
            float[] output = new float[genes.size() / 4];
            for (int i = 0; i < output.length; i++) {
                output[i] = bytebuffer.getFloat();
            }
            ln += "\t" + Strings.concat(output, Strings.tab);

            for (int g = 0; g < genes.size(); g++) {
                ln2 += "\t" + snp.n;
            }
            tfout.writeln(ln);
            tfout2.writeln(ln2);
            pb.iterate();
        }
        pb.close();
        bf.close();
        tfout.close();
        tfout2.close();
    }


}
