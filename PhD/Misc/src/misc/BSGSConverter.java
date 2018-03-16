/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class BSGSConverter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            ProbeTranslation pb = new ProbeTranslation();
            HashMap<String, String> ht12ToMetaId = pb.getProbeTranslation("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", "HumanHT-12_V4_0_R1_15002873_B.txt", "Probe");

            for (int i = 0; i < 11; i++) {
                int pvalCol = 2;
                int zsCol = 3;
                String fName = "eQTLs.txt";

                if (i > 0) {
                    fName = "PermutedEQTLsPermutationRound" + i + ".txt.gz";
                    pvalCol = i + 8;
                    zsCol = i + 18;
                }

                TextFile in = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/Replication/BSGS/BSGS_trans_eQTL_replication_results2.csv", TextFile.R);
                TextFile out = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/Replication/BSGS/ParsedUnsorted/" + fName, TextFile.W);
                String header = in.readLine();


                // PValue	SNPName	SNPChr	SNPChrPos	ProbeName	ProbeChr	ProbeCenterChrPos	CisTrans	SNPType	AlleleAssessed	OverallZScore	DatasetsWhereSNPProbePairIsAvailableAndPassesQC	DatasetsZScores	DatasetsNrSamples	IncludedDatasetsMeanProbeExpression	IncludedDatasetsProbeExpressionVariance	HGNCName	IncludedDatasetsCorrelationCoefficient	Meta-Beta (SE)	Beta (SE)	FoldChange	FDR

                out.writeln("PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tIncludedDatasetsMeanProbeExpression\tIncludedDatasetsProbeExpressionVariance\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta (SE)\tBeta (SE)\tFoldChange\tFDR");

                // SNPName,ProbeName,PValue,ZScore,ZScoreDirAllele,MinorAllele,SNPType,SampleN,MAF,Permutation1Pval,Permutation2Pval,Permutation3Pval,Permutation4Pval,Permutation5Pval,Permutation6Pval,Permutation7Pval,Permutation8Pval,Permutation9Pval,Permutation10Pval,Permutation1ZScore,Permutation2ZScore,Permutation3ZScore,Permutation4ZScore,Permutation5ZScore,Permutation6ZScore,Permutation7ZScore,Permutation8ZScore,Permutation9ZScore,Permutation10ZScore
                String[] elems = in.readLineElems(TextFile.comma);
                while (elems != null) {

                    String[] data = new String[22];

                    data[0] = elems[pvalCol];
                    data[1] = elems[0];
                    data[4] = ht12ToMetaId.get(elems[1]);
                    data[8] = elems[6];
                    data[9] = elems[4];
                    data[10] = elems[zsCol];

                    String output = Strings.concat(data, Strings.tab);
                    out.writeln(output);
                    elems = in.readLineElems(TextFile.comma);
                }

                out.close();
                in.close();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
