package nl.harmjanwestra.playground.biogen.mrvisualisation;

import com.itextpdf.text.DocumentException;
import umcg.genetica.graphics.Grid;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.*;

public class MRViz {

    public static void main(String[] args) {

        String mrfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-01-27-MRViz\\AllMRData.txt";
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-01-27-MRViz\\AllMRData-OR.pdf";


        MRViz MRViz = new MRViz();
        try {
            boolean makeEQTLPositive = true;
            mrfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-01-27-MRViz\\AllMRData-DiseaseRisk.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-01-27-MRViz\\AllMRData-OR.pdf";
//            MRViz.runOR(mrfile, output, true, makeEQTLPositive);
//
            mrfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-01-27-MRViz\\AllMRData-Quantitative.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-01-27-MRViz\\AllMRData-Beta.pdf";
//            MRViz.runOR(mrfile, output, false, makeEQTLPositive);


//            mrfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2022-05-MRViz\\table1-ms.txt";
//            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2022-05-MRViz\\table1-ms.pdf";
//            MRViz.runOR(mrfile, output, true, makeEQTLPositive);

            mrfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2022-05-MRViz\\suptable12-all-or.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2022-05-MRViz\\suptable12-all-or.pdf";
            MRViz.runOR(mrfile, output, true, makeEQTLPositive);

            mrfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2022-05-MRViz\\suptable12-all-quant.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2022-05-MRViz\\suptable12-all-quant.pdf";
            MRViz.runOR(mrfile, output, false, makeEQTLPositive);

        } catch (IOException e) {
            e.printStackTrace();
        } catch (DocumentException e) {
            e.printStackTrace();
        }

    }

    public void run(String mrFile, String output) throws IOException, DocumentException {

        TextFile tf = new TextFile(mrFile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        HashMap<String, ArrayList<MREvent>> mrEvents = new HashMap<>();

        while (elems != null) {
            // SNP	EA	NONEA	gene	beta_exposure	se_exposure	outcome	beta_outcome	se_outcome	WR	se	default.coloc.PP4	nsnps.coloc.PP4
            if (elems.length > 12) {
                MREvent event = new MREvent();
                event.snp = elems[0];
                event.EA = elems[1];
                event.NONEA = elems[2];
                event.gene = elems[3];
                event.betaExposure = Double.parseDouble(elems[4]);
                event.seExposure = Double.parseDouble(elems[5]);
                event.outcome = elems[6];
                event.betaOutcome = Double.parseDouble(elems[7]);
                event.seOutcome = Double.parseDouble(elems[8]);
                event.waldratio = Double.parseDouble(elems[9]);
                event.seWaldratio = Double.parseDouble(elems[10]);
                event.pp4def = Double.parseDouble(elems[11]);
                event.pp4bon = Double.parseDouble(elems[12]);
                ArrayList<MREvent> outcomeEvents = mrEvents.get(event.outcome);
                if (outcomeEvents == null) {
                    outcomeEvents = new ArrayList<>();
                }
                outcomeEvents.add(event);
                mrEvents.put(event.outcome, outcomeEvents);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        for (Map.Entry<String, ArrayList<MREvent>> pair : mrEvents.entrySet()) {
            Collections.sort(pair.getValue());
        }

        Grid grid = new Grid(Grid.SIZE.LETTER, 1000, 1500, 1, 1, 100, 100);
        System.out.println("plotting");
        MRPanel panel = new MRPanel(1, 1);
        panel.determineRangeOverAllEvents();
        panel.overlapBoxes();
        grid.addPanel(panel);
        panel.setData(mrEvents);

        grid.draw(output);
    }

    public void runOR(String mrFile, String output, boolean applylog, boolean makeEQTLPositive) throws IOException, DocumentException {

        TextFile tf = new TextFile(mrFile, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        int snpcol = -1;
        int eacol = -1;
        int noneacol = -1;
        int genecol = -1;
        int betaExpCol = -1;
        int seExpCol = -1;
        int outcomeCol = -1;
        int betaOutCol = -1;
        int seOutCol = -1;
        int wrcol = -1;
        int sewrcol = -1;
        int pwrcol = -1;
//        int pp4defcol = -1;
//        int pp4boncol = -1;

        // chrom	pos	SNP	EA	NONEA	proxy_SNP	gene	beta_exposure	se_exposure	p_exposure	outcome	beta_outcome	se_outcome	p_outcome	WR	se	p
        // id	chrom	pos	SNP	EA	NONEA	proxy_used	proxy_SNP	ENSG	gene	beta_exposure	se_exposure	p_exposure	outcome	MRBase_id	beta_outcome	se_outcome	p_outcome	WR	se	p	pass bonferroni correction	default.coloc.PP4	nsnps.coloc.PP4	coloc.relaxed	coloc.strict	Astrocyte_Beta	EndothelialCell_Beta	Macrophage_Beta	Neuron_Beta	Oligodendrocyte_Beta	Astrocyte_FDR	EndothelialCell_FDR	Macrophage_FDR	Neuron_FDR	Oligodendrocyte_FDR	deconQTLFlip	n-signif.	DDG2P (accessed: 20-12-16)	DDG2P - mutation consequence	OrphaNet

        // SNP	Gene	Beta Exposure	SE Exposure	Outcome	Beta Outcome	SE Outcome	WaldRatio	WaldRatio SE	COLOC PP4 Corrected for number of SNPs
        // ID	Chromosome	Position	SNP	Effect Allele	Non Effect Allele	Proxy Used	Proxy SNP	Ensembl Gene ID	Gene	Beta Exposure	SE Exposure	P Exposure	Outcome	MRBase ID
        // Beta Outcome	SE Outcome	P Outcome	WaldRatio	WaldRatio SE	WaldRatio P	Passes Bonferroni correction	Default Coloc PP4
        // COLOC PP4 Corrected for number of SNPs	Coloc Relaxed	Coloc Strict
        //
        // Ensembl Gene ID	Gene	Beta Exposure	SE Exposure	P Exposure	Outcome	MRBase ID	Beta Outcome	SE Outcome	P Outcome	WaldRatio	WaldRatio SE	WaldRatio P	Passes Bonferroni correction	Default Coloc PP4	COLOC PP4 Corrected for number of SNPs	Coloc Relaxed	Coloc Strict
        for (int i = 0; i < header.length; i++) {
            String q = header[i];
            if (q.equals("SNP")) {
                snpcol = i;
            } else if (q.equals("EA") || q.equals("Effect Allele")) {
                eacol = i;
            } else if (q.equals("NONEA") || q.equals("Non Effect Allele")) {
                noneacol = i;
            } else if (q.equals("gene") || q.equals("Gene")) {
                genecol = i;
            } else if (q.equals("beta_exposure") || q.equals("Beta Exposure")) {
                betaExpCol = i;
            } else if (q.equals("se_exposure") || q.equals("SE Exposure")) {
                seExpCol = i;
            } else if (q.equals("outcome") || q.equals("Outcome")) {
                outcomeCol = i;
            } else if (q.equals("beta_outcome") || q.equals("Beta Outcome")) {
                betaOutCol = i;
            } else if (q.equals("se_outcome") || q.equals("SE Outcome")) {
                seOutCol = i;
            } else if (q.equals("WR") || q.equals("WaldRatio")) {
                wrcol = i;
            } else if (q.equals("se") || q.equals("WaldRatio SE")) {
                sewrcol = i;
            } else if (q.equals("p") || q.equals("WaldRatio P")) {
                pwrcol = i;
            }

//            else if (q.equals("")) {
//
//            } else if (q.equals("")) {
//
//            } else if (q.equals("")) {
//
//            }
        }

        String[] elems = tf.readLineElems(TextFile.tab);
        HashMap<String, ArrayList<MREvent>> mrEvents = new HashMap<>();

        while (elems != null) {
            // SNP	EA	NONEA	gene	beta_exposure	se_exposure	outcome	beta_outcome	se_outcome	WR	se	default.coloc.PP4	nsnps.coloc.PP4
            if (elems.length > 12) {
                MREvent event = new MREvent();
                event.snp = elems[snpcol];
                event.EA = elems[eacol];
                event.NONEA = elems[noneacol];
                event.gene = elems[genecol];
                event.betaExposure = Double.parseDouble(elems[betaExpCol]);
                event.seExposure = Double.parseDouble(elems[seExpCol]);
                event.outcome = elems[outcomeCol];
                event.betaOutcome = Double.parseDouble(elems[betaOutCol]);
                event.seOutcome = Double.parseDouble(elems[seOutCol]);
                event.pWR = elems[pwrcol];
                if (makeEQTLPositive) {
                    if (event.betaExposure < 0) {
                        event.betaExposure *= -1;
                        event.betaOutcome *= -1;
                        event.EA = elems[eacol];
                        event.NONEA = elems[noneacol];
                    }
                }

                event.waldratio = Double.parseDouble(elems[wrcol]);
                event.seWaldratio = Double.parseDouble(elems[sewrcol]);
//                event.pp4def = Double.parseDouble(elems[pp4defcol]);
//                event.pp4bon = Double.parseDouble(elems[pp4boncol]);
                ArrayList<MREvent> outcomeEvents = mrEvents.get(event.outcome);
                if (outcomeEvents == null) {
                    outcomeEvents = new ArrayList<>();
                }
                outcomeEvents.add(event);
                mrEvents.put(event.outcome, outcomeEvents);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        for (Map.Entry<String, ArrayList<MREvent>> pair : mrEvents.entrySet()) {
            Collections.sort(pair.getValue());
        }

        Grid grid = new Grid(Grid.SIZE.LETTER, 1000, 1500, 1, 1, 100, 100);
        System.out.println("plotting");
        MRPanelOR panel = new MRPanelOR(1, 1);
        panel.applyLogToOutcome(applylog);
        // panel.determineRangeOverAllEvents();
        // panel.overlapBoxes();
        grid.addPanel(panel);
        panel.setData(mrEvents);

        grid.draw(output);
    }

}
