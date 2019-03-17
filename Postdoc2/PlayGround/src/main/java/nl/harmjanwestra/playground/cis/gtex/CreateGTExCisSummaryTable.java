package nl.harmjanwestra.playground.cis.gtex;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

public class CreateGTExCisSummaryTable {

    public static void main(String[] args) {

        String infile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-06-05-cis-gtexrepl\\usingFullGTExResults\\cisrepltable.txt";

        String trans = "D:\\Sync\\GDrive\\FrankeSwertzLab\\Projects\\eQTLgen\\DataFreezes\\2018-07-26-transEQTL-Replication-MasterTable.txt.gz";

        CreateGTExCisSummaryTable c = new CreateGTExCisSummaryTable();

        try {
            // c.gtfToProbeAnnotationFile(infile, false);
            c.run(trans, true);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String infile, boolean useFDRCol) throws IOException {


        TextFile tf = new TextFile(infile, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        ArrayList<String> tissues = new ArrayList<>();
        for (String s : header) {
            if (!s.equals("-") && !s.equals("eQTLGen") && s.length() > 0) {
                tissues.add(s);
            }
        }
        String[] header2 = tf.readLineElems(TextFile.tab);

        ArrayList<Integer> zcols = new ArrayList<>();
        ArrayList<Integer> sigcols = new ArrayList<>();
        for (int i = 0; i < header2.length; i++) {
            if (header2[i].equals("Z")) {
                zcols.add(i);
            }
            if (!useFDRCol) {
                if (header2[i].equals("Significant")) {
                    sigcols.add(i);
                }
            } else {
                if (header2[i].equals("FDR")) {
                    sigcols.add(i);
                }
            }

        }

        String[] elems = tf.readLineElems(TextFile.tab);

        int[] shared = new int[zcols.size()];
        int[] sharedsamedir = new int[zcols.size()];
        int[] sharedsig = new int[zcols.size()];
        int[] sharedsigsamedir = new int[zcols.size()];

        int total = 0;
        while (elems != null) {

            Double ref = null;
            boolean refSig = false;
            for (int i = 0; i < zcols.size(); i++) {
                Integer colid = zcols.get(i);
                Integer sigcol = sigcols.get(i);
                if (i == 0) {
                    ref = Double.parseDouble(elems[colid]);
                    refSig = Boolean.parseBoolean(elems[sigcol]);
                } else {
                    if (!elems[colid].equals("-")) {
                        double z = Double.parseDouble(elems[colid]);
                        boolean zsig = false;
                        if (!useFDRCol) {
                            zsig = Boolean.parseBoolean(elems[sigcol]);
                        } else {
                            Double fdr = Double.parseDouble(elems[sigcol]);
                            if (fdr < 0.05) {
                                zsig = true;
                            }
                        }
                        shared[i]++;
                        if (zsig) {
                            sharedsig[i]++;
                        }
                        if (z * ref > -1) {
                            sharedsamedir[i]++;
                            if (zsig) {
                                sharedsigsamedir[i]++;
                            }
                        }

                    }

                }
            }
            total++;
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        for (int i = 0; i < tissues.size(); i++) {
            System.out.println(tissues.get(i) + "\t" + shared[i + 1] + "\t" + sharedsamedir[i + 1] + "\t" + sharedsig[i + 1] + "\t" + sharedsigsamedir[i + 1]);
        }
        System.out.println(total);
    }
}
