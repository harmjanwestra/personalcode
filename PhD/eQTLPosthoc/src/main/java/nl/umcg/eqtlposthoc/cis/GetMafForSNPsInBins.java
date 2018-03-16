/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.cis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

/**
 *
 * @author harmjan
 */
public class GetMafForSNPsInBins {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            GetMafForSNPsInBins s = new GetMafForSNPsInBins();
            s.run();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    private TriTyperGenotypeData ds;
    private SNPLoader loader;
    private TextFile outfile;
    int snpcounter = 0;
    private HashMap<Integer, ArrayList<Double>> frequenciesPerBin;
    private HashMap<Integer, ArrayList<Double>> frequenciesPerBinRealData;

    public void run() throws IOException {

        ds = new TriTyperGenotypeData();
        ds.load("/Data/GeneticalGenomicsDatasets/BloodHT12ImputeTriTyper/");
        loader = ds.createSNPLoader();



        frequenciesPerBin = new HashMap<Integer, ArrayList<Double>>();
        frequenciesPerBinRealData = new HashMap<Integer, ArrayList<Double>>();

        String indir = "/Volumes/iSnackHD/Skydrive/Cis-SNPs/";

        for (int perm = 0; perm < 10; perm++) {
            for (int iter = 0; iter < 16000; iter += 1000) {
                String iterStr = "Bin" + iter + "-" + (iter + 1000) + ".txt";
                String binfile = "";
                String proxyfile = "";
                if (perm == 0) {
                    binfile = indir + "RealData/" + iterStr;
                    proxyfile = "/Volumes/iSnackHD/Skydrive/Cis-SNPs/RealData/AllSNPs.txt-WithProxies.txt";
                } else {
                    proxyfile = "/Volumes/iSnackHD/Skydrive/Cis-SNPs/AllSNPsPerm.txt-WithProxies.txt";
                    binfile = indir + "PermutationRound" + perm + "/" + iterStr;
                }
                this.getMaf(perm, iter, proxyfile, binfile);
//                System.out.println(perm + "\t" + iter + "\t" + snpcounter);
            }
        }

        int maxnr = 0;


        int[][] bins = new int[26][16];
        int[] bintotals = new int[16];

        for (int iter = 0; iter < 16000; iter += 1000) {

            ArrayList<Double> freqPerm = frequenciesPerBinRealData.get(iter);
            if (freqPerm.size() > maxnr) {
                maxnr = freqPerm.size();
            }

            int pnum = iter / 1000;
            bintotals[pnum] = freqPerm.size();
            for (Double d : freqPerm) {
                int binNo = (int) Math.floor(d * 50d);

                if (binNo > 50) {
                    System.out.println(d);
                }
                bins[binNo][pnum]++;
            }
        }

        outfile = new TextFile("/Volumes/iSnackHD/Skydrive/Cis-SNPs/MafForSNPsRealFreqNrs.txt", TextFile.W);
        String out = "MAF";
        for (int i = 0; i < 16000; i += 1000) {
            out += "\t" + (i + 1000);
        }
        outfile.writeln(out);
        for (int i = 0; i < bins.length; i++) {
            out = ((double) i * 2 / 100) + "";
            for (int j = 0; j < bins[i].length; j++) {
//                out += "\t" + ((double) bins[i][j] / bintotals[j]);
                out += "\t" + bins[i][j];
            }
            outfile.writeln(out);
        }
        outfile.close();

//        for (int i = 0; i < maxnr; i++) {
//            String output = "" + i;
//            for (int iter = 0; iter < 16000; iter += 1000) {
//                ArrayList<Double> freqPerm = frequenciesPerBinRealData.get(iter);
//                if (i > freqPerm.size() - 1) {
//                    output += "\t";
//                } else {
//                    output += "\t" + freqPerm.get(i);
//                }
//            }
//            outfile.writeln(output);
//        }
//        outfile.close();


        maxnr = 0;
        bins = new int[26][16];
        bintotals = new int[16];
        for (int iter = 0; iter < 16000; iter += 1000) {
            ArrayList<Double> freqPerm = frequenciesPerBin.get(iter);
            if (freqPerm.size() > maxnr) {
                maxnr = freqPerm.size();
            }
            int pnum = iter / 1000;
            bintotals[pnum] = freqPerm.size();
            for (Double d : freqPerm) {
                int binNo = (int) Math.floor(d * 50d);
                if (binNo > 50) {
                    System.out.println(d);
                }
                bins[binNo][pnum]++;
            }
        }


        outfile = new TextFile("/Volumes/iSnackHD/Skydrive/Cis-SNPs/MafForSNPsFreqNrs.txt", TextFile.W);
        out = "MAF";
        for (int i = 0; i < 16000; i += 1000) {
            out += "\t" + (i + 1000);
        }
        outfile.writeln(out);
        for (int i = 0; i < bins.length; i++) {
            out = i + "";
            for (int j = 0; j < bins[i].length; j++) {
//                out += "\t" + ((double) bins[i][j] / bintotals[j]);
                out += "\t" + bins[i][j];
            }
            outfile.writeln(out);
        }
        outfile.close();

    }

    public void getMaf(int perm, int binno, String proxyfile, String binfile) throws IOException {
        ArrayList<Double> freq = null;
        if (perm == 0) {
            freq = frequenciesPerBinRealData.get(binno);
        } else {
            freq = frequenciesPerBin.get(binno);
        }
        if (freq == null) {
            freq = new ArrayList<Double>();
        }

        HashSet<String> snpsToQuery = new HashSet<String>();
        TextFile in = new TextFile(binfile, TextFile.R);
        snpsToQuery.addAll(in.readAsArrayList());
        in.close();

        TextFile in2 = new TextFile(proxyfile, TextFile.R);
        String[] elems = in2.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[0];
            String proxy = elems[1];
            if (snpsToQuery.contains(snp)) {
                snpsToQuery.add(proxy);
            }
            elems = in2.readLineElems(TextFile.tab);
        }

        in.close();

        for (String snp : snpsToQuery) {
            Integer snpId = ds.getSnpToSNPId().get(snp);
            if (snpId != null) {
                snpcounter++;
                SNP snpObj = ds.getSNPObject(snpId);
                loader.loadGenotypes(snpObj);
                double maf = snpObj.getMAF();
//                System.out.println(perm + "\t" + binno + "\t" + snp + "\t" + maf);
                freq.add(maf);
                snpObj.clearGenotypes();
            }
        }

        if (perm == 0) {
            frequenciesPerBinRealData.put(binno, freq);
        } else {
            frequenciesPerBin.put(binno, freq);
        }
    }
}
