/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class eQTLsForTraits {

    HashSet<String> allowedTraits = null;
    HashMap<String, String> snpToChr = new HashMap<String, String>();
    HashMap<String, String> snpToPos = new HashMap<String, String>();

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        eQTLsForTraits e = new eQTLsForTraits();
        e.run();
    }
    private GWASCatalog catalog;

    private void loadDBSNP(String dbSNP) throws IOException {

        ArrayList<GWASSNP> gwassnps = new ArrayList<GWASSNP>();
        gwassnps.addAll(catalog.getSnps());

        HashSet<String> gwasSNPstr = new HashSet<String>();
        for (int i = 0; i < gwassnps.size(); i++) {
            gwasSNPstr.add(gwassnps.get(i).getName());
        }

        TextFile dbsnptf = new TextFile(dbSNP, TextFile.R);
        String[] dbsnpelems = dbsnptf.readLineElems(TextFile.tab);

        while (dbsnpelems != null) {

            String snpstr = dbsnpelems[2];
            if (gwasSNPstr.contains(snpstr)) {
                snpToChr.put(snpstr, dbsnpelems[0]);
                snpToPos.put(snpstr, dbsnpelems[1]);
            }
            dbsnpelems = dbsnptf.readLineElems(TextFile.tab);
        }

        dbsnptf.close();
    }

    public void run() {
        try {


            catalog = new GWASCatalog();

            catalog.read("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog.txt");

//	    loadTraitSelection("");

//	    GWASTrait[] selectedTraitsFromCatalog = c.getTraitsForCertainKey("body mass index");


            HashSet<GWASSNP> selectedSNPsFromCatalog = new HashSet<GWASSNP>();

//	    selectedSNPsFromCatalog.addAll(Arrays.asList(catalog.getSNPsForTraitContainingKey("Diabetes related insulin traits")));
//	    selectedSNPsFromCatalog.addAll(Arrays.asList(catalog.getSNPsForTraitContainingKey("Diabetes (incident)")));
//	    selectedSNPsFromCatalog.addAll(Arrays.asList(catalog.getSNPsForTraitContainingKey("Type 2 diabetes")));
//	    selectedSNPsFromCatalog.addAll(Arrays.asList(catalog.getSNPsForTraitContainingKey("Type 2 diabetes and other traits")));

            selectedSNPsFromCatalog.addAll(Arrays.asList(catalog.getSnps().toArray(new GWASSNP[0])));

            System.out.println(selectedSNPsFromCatalog.size());



            loadDBSNP("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-27-SNPMappings-dbSNP130.txt.gz");







            HashMap<GWASSNP, HashSet<String>> cisEQTLsForSNPs = new HashMap<GWASSNP, HashSet<String>>();
            HashMap<GWASSNP, HashSet<String>> transEQTLsForSNPs = new HashMap<GWASSNP, HashSet<String>>();

            TextFile eqtlfile = new TextFile("/Volumes/iSnackHD/Data/Projects/LudeFranke/eqtlresults/2012-06-14-CISTRANS-Reactome_BPMean-DiseaseSNPs-CisRegressedOut-META/eQTLProbesFDR0.05.txt", TextFile.R);
            String[] elems = eqtlfile.readLineElems(TextFile.tab);
            while (elems != null) {

                String snp = elems[1];
                String gene = elems[4];

                GWASSNP s = catalog.getSnpToObj().get(snp);
                if (s != null) {
                    if (!gene.equals("-")) {
                        if (selectedSNPsFromCatalog.contains(s)) {

                            boolean ciseffect = false;

                            if (elems[2].equals(elems[5])) {
                                Integer snppos = Integer.parseInt(elems[3]);
                                Integer probepos = Integer.parseInt(elems[6]);

                                if (Math.abs(snppos - probepos) < 5000000) {
                                    ciseffect = true;
                                }

                            }

                            if (ciseffect) {
                                HashSet<String> cis = cisEQTLsForSNPs.get(s);
                                if (cis == null) {
                                    cis = new HashSet<String>();
                                }

                                cis.add(gene);

                                cisEQTLsForSNPs.put(s, cis);
                            } else {
                                HashSet<String> trans = transEQTLsForSNPs.get(s);

                                if (trans == null) {
                                    trans = new HashSet<String>();
                                }

                                trans.add(gene);
                                transEQTLsForSNPs.put(s, trans);
                            }




                        }
                    }

                }


                elems = eqtlfile.readLineElems(TextFile.tab);
            }


            eqtlfile.close();


            eqtlfile = new TextFile("/Volumes/iSnackHD/Data/Projects/LudeFranke/eqtlresults/2012-06-14-CISTRANS-Reactome_BPMean-DiseaseSNPs-CisRegressedOut-META/eQTLProbesFDR0.05.txt", TextFile.R);
            elems = eqtlfile.readLineElems(TextFile.tab);
            while (elems != null) {

                String snp = elems[1];
                String gene = elems[4];

                GWASSNP s = catalog.getSnpToObj().get(snp);

                if (s != null) {
                    if (!gene.equals("-")) {
                        if (selectedSNPsFromCatalog.contains(s)) {

                            boolean ciseffect = false;

                            if (elems[2].equals(elems[5])) {
                                Integer snppos = Integer.parseInt(elems[3]);
                                Integer probepos = Integer.parseInt(elems[6]);

                                if (Math.abs(snppos - probepos) < 5000000) {
                                    ciseffect = true;
                                }

                            }

                            if (ciseffect) {
                                HashSet<String> cis = cisEQTLsForSNPs.get(s);
                                if (cis == null) {
                                    cis = new HashSet<String>();
                                }

                                cis.add(gene);

                                cisEQTLsForSNPs.put(s, cis);
                            } else {
                                HashSet<String> trans = transEQTLsForSNPs.get(s);

                                if (trans == null) {
                                    trans = new HashSet<String>();
                                }

                                trans.add(gene);
                                transEQTLsForSNPs.put(s, trans);
                            }




                        }
                    }
                }






                elems = eqtlfile.readLineElems(TextFile.tab);
            }


            eqtlfile.close();


            GWASSNP[] snparr = selectedSNPsFromCatalog.toArray(new GWASSNP[selectedSNPsFromCatalog.size()]);
            TextFile out = new TextFile("/Volumes/iSnackHD/SkyDrive/PathwaysInfluencedByTraitAssocSNPs.txt", TextFile.W);
            for (GWASSNP s : snparr) {

                HashSet<String> cis = cisEQTLsForSNPs.get(s);
                HashSet<String> trans = transEQTLsForSNPs.get(s);

                if (s.getChr() != 6) {
                    if (cis != null && trans != null) {
                        GWASTrait[] traitsAssoc = s.getAssociatedTraitsArray();
                        String[] str = new String[traitsAssoc.length];
                        for (int q = 0; q < str.length; q++) {
                            str[q] = traitsAssoc[q].getName();
                        }

                        String traitStr = Strings.concat(str, Strings.comma);

                        String[] cisFX = cis.toArray(new String[0]);
                        String[] transFX = trans.toArray(new String[0]);

                        String cisConcat = Strings.concat(cisFX, Strings.comma);
                        String transConcat = Strings.concat(transFX, Strings.comma);

                        System.out.println(traitStr + "\t" + s.getName() + "\t" + snpToChr.get(s.getName()) + "\t" + cisConcat + "\t" + transConcat);

                        out.writeln(traitStr + "\t" + s.getName() + "\t" + snpToChr.get(s.getName()) + "\t" + cisFX.length + "\t" + transFX.length + "\t" + cisConcat + "\t" + transConcat);
                    }
                }

            }
            out.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void loadTraitSelection(String traitselectionfile) throws IOException {
        System.out.println("Loading trait selection from: " + traitselectionfile);
        TextFile tf = new TextFile(traitselectionfile, TextFile.R);
        allowedTraits = new HashSet<String>();
        String[] elems = tf.readLineElems(TextFile.tab);

        while (elems != null) {
            if (elems.length > 1) {
                if (elems[1].equals("FALSE")) {
                    allowedTraits.add(elems[0]);
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        System.out.println(allowedTraits.size() + " traits selected");
    }
}
