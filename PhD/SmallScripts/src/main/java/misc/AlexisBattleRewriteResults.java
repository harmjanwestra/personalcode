/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author harmjan
 */
public class AlexisBattleRewriteResults {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            TextFile tfE = new TextFile("/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.5-HT12v3.txt", TextFile.R);

            tfE.readLine();
            String[] tfElems = tfE.readLineElems(TextFile.tab);
            HashMap<String, String> comboToGene = new HashMap<String, String>();
            while (tfElems != null) {
                String[] genelems = tfElems[16].split(",");
                for (String gene : genelems) {
                    comboToGene.put(tfElems[1] + "-" + gene, tfElems[4]);
                }

                System.out.println(tfElems[1] + "-" + tfElems[16] + "\t" + tfElems[4]);
                tfElems = tfE.readLineElems(TextFile.tab);
            }
            tfE.close();

            TextFile tf = new TextFile("/Volumes/iSnackHD/Data/Projects/AlexisBattle/2013-11-07-TransRepl/westra_trans_replication.txt", TextFile.R);
            TextFile tfout = new TextFile("/Volumes/iSnackHD/Data/Projects/AlexisBattle/2013-11-07-TransRepl/westra_trans_replication-eQTLFile-DirectionInverted.txt", TextFile.W);

            Correlation.correlationToZScore(922);

            System.out.println("");
            System.out.println("");
            System.out.println("");

            tf.readLine();
            tfout.writeln("PValue	SNPName	SNPChr	SNPChrPos	ProbeName	ProbeChr	ProbeCenterChrPos	CisTrans	SNPType	AlleleAssessed	OverallZScore	DatasetsWhereSNPProbePairIsAvailableAndPassesQC	DatasetsZScores	DatasetsNrSamples	IncludedDatasetsMeanProbeExpression	IncludedDatasetsProbeExpressionVariance	HGNCName	IncludedDatasetsCorrelationCoefficient	Meta-Beta (SE)	Beta (SE)	FoldChange	FDR");
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                double correlation = -Double.parseDouble(elems[4]);
                Double z = Correlation.convertCorrelationToZScore(922, correlation);
                double p = ZScores.zToP(z);
                String probe = comboToGene.get(elems[0] + "-" + elems[1]);
                System.out.println(elems[0] + "-" + elems[1] + "\t" + probe);
                tfout.writeln(p + "\t" + elems[0] + "\t0\t0\t" + probe + "\t0\t0\ttrans\t" + elems[2].replaceAll(",", "/") + "\t" + elems[3] + "\t" + z + "\t0\t" + z + "\t922\t0\t0\t" + elems[1] + "\t" + elems[4] + "\t0\t0\t0\t" + elems[7]);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            tfout.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /*
     PValue	SNPName	SNPChr	SNPChrPos	ProbeName	ProbeChr	ProbeCenterChrPos	CisTrans	SNPType	AlleleAssessed	OverallZScore	DatasetsWhereSNPProbePairIsAvailableAndPassesQC	DatasetsZScores	DatasetsNrSamples	IncludedDatasetsMeanProbeExpression	IncludedDatasetsProbeExpressionVariance	HGNCName	IncludedDatasetsCorrelationCoefficient	Meta-Beta (SE)	Beta (SE)	FoldChange	FDR	FDR
     2.0013718859520384E-16	rs2158354	21	42710175	6270743	21	42718970	cis	T/C	C	8.222082216130435	BloodHT12v3-WithCovariate	8.222082216130435	1240	1.7548686518268602E-17	0.3692268414037504	FAM3B	0.11958462761054467	0.0 (0.0)	0.0 (0.0)	0.0	0.0
     */
}
