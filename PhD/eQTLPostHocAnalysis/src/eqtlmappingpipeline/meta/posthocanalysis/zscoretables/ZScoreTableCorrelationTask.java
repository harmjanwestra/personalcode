/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.zscoretables;

import java.util.concurrent.Callable;
import java.util.concurrent.LinkedBlockingQueue;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.util.RankDoubleArray;

/**
 *
 * @author harmjan
 */
public class ZScoreTableCorrelationTask implements Callable<Pair<Integer, Double>> {

    DoubleMatrixDataset<String, String> pathwayZScoreData;
    DoubleMatrixDataset<String, String> eQTLZScoreData;
    int sharedgenes;
    int[] colToCol;
    int correlationdistributionr2length;
    int i;
    int j;
    private final boolean nonparametric;
    RankDoubleArray rda = new RankDoubleArray();

    ZScoreTableCorrelationTask(DoubleMatrixDataset<String, String> pathwayZScoreData, int sharedgenes, DoubleMatrixDataset<String, String> eQTLZScoreData, int[] colToCol, int i, int j, int length, boolean nonparametric) {
	this.pathwayZScoreData = pathwayZScoreData;
	this.sharedgenes = sharedgenes;
	this.eQTLZScoreData = eQTLZScoreData;
	this.colToCol = colToCol;
	this.i = i;
	this.j = j;
	this.correlationdistributionr2length = length;

	this.nonparametric = nonparametric;

    }

    @Override
    public Pair<Integer, Double> call() throws Exception {


//	for (int j = 0; j < pathwayZScoreData.nrRows; j++) {
	double[] valsX = new double[sharedgenes];
	double[] valsY = new double[sharedgenes];
	int itr = 0;
	for (int p = 0; p < eQTLZScoreData.nrCols; p++) {
	    if (colToCol[p] != -1) {
		valsX[itr] = Math.abs(eQTLZScoreData.rawData[i][p]);
		valsY[itr] = pathwayZScoreData.rawData[j][colToCol[p]];
		itr++;
	    }
	}

	if (nonparametric) {

	    valsX = rda.rank(valsX);
	    valsY = rda.rank(valsY);
	}

	double correlation = JSci.maths.ArrayMath.correlation(valsX, valsY);

	double r2 = correlation * correlation;


	return new Pair<Integer, Double>(j, r2);

//	}
    }
}
