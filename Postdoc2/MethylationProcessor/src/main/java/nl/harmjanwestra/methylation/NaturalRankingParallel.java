package nl.harmjanwestra.methylation;

import org.apache.commons.math3.exception.MathInternalError;
import org.apache.commons.math3.exception.NotANumberException;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.IntStream;

public class NaturalRankingParallel implements RankingAlgorithm {
	public static final NaNStrategy DEFAULT_NAN_STRATEGY;
	public static final TiesStrategy DEFAULT_TIES_STRATEGY;
	private final NaNStrategy nanStrategy;
	private final TiesStrategy tiesStrategy;
	private final RandomDataGenerator randomData;

	public NaturalRankingParallel() {
		this.tiesStrategy = DEFAULT_TIES_STRATEGY;
		this.nanStrategy = DEFAULT_NAN_STRATEGY;
		this.randomData = null;
	}

	public NaturalRankingParallel(TiesStrategy tiesStrategy) {
		this.tiesStrategy = tiesStrategy;
		this.nanStrategy = DEFAULT_NAN_STRATEGY;
		this.randomData = new RandomDataGenerator();
	}

	public NaturalRankingParallel(NaNStrategy nanStrategy) {
		this.nanStrategy = nanStrategy;
		this.tiesStrategy = DEFAULT_TIES_STRATEGY;
		this.randomData = null;
	}

	public NaturalRankingParallel(NaNStrategy nanStrategy, TiesStrategy tiesStrategy) {
		this.nanStrategy = nanStrategy;
		this.tiesStrategy = tiesStrategy;
		this.randomData = new RandomDataGenerator();
	}

	public NaturalRankingParallel(RandomGenerator randomGenerator) {
		this.tiesStrategy = TiesStrategy.RANDOM;
		this.nanStrategy = DEFAULT_NAN_STRATEGY;
		this.randomData = new RandomDataGenerator(randomGenerator);
	}

	public NaturalRankingParallel(NaNStrategy nanStrategy, RandomGenerator randomGenerator) {
		this.nanStrategy = nanStrategy;
		this.tiesStrategy = TiesStrategy.RANDOM;
		this.randomData = new RandomDataGenerator(randomGenerator);
	}

	public NaNStrategy getNanStrategy() {
		return this.nanStrategy;
	}

	public TiesStrategy getTiesStrategy() {
		return this.tiesStrategy;
	}

	public double[] rank(double[] data) {
		IntDoublePair[] ranks = new IntDoublePair[data.length];

		IntDoublePair[] finalRanks = ranks;
		IntStream.range(0, data.length).parallel().forEach(i -> {
			finalRanks[i] = new IntDoublePair(data[i], i);
		});
		ranks = finalRanks;

		List<Integer> nanPositions = null;
		switch (this.nanStrategy) {
			case MAXIMAL:
				this.recodeNaNs(ranks, 1.0D / 0.0);
				break;
			case MINIMAL:
				this.recodeNaNs(ranks, -1.0D / 0.0);
				break;
			case REMOVED:
				ranks = this.removeNaNs(ranks);
				break;
			case FIXED:
				nanPositions = this.getNanPositions(ranks);
				break;
			case FAILED:
				nanPositions = this.getNanPositions(ranks);
				if (nanPositions.size() > 0) {
					throw new NotANumberException();
				}
				break;
			default:
				throw new MathInternalError();
		}

		Arrays.parallelSort(ranks);

		double[] out = new double[ranks.length];
		int pos = 1;
		out[ranks[0].getPosition()] = (double) pos;
		List<Integer> tiesTrace = new ArrayList();
		tiesTrace.add(ranks[0].getPosition());

		for (int i = 1; i < ranks.length; ++i) {
			if (Double.compare(ranks[i].getValue(), ranks[i - 1].getValue()) > 0) {
				pos = i + 1;
				if (tiesTrace.size() > 1) {
					this.resolveTie(out, tiesTrace);
				}

				tiesTrace = new ArrayList();
				tiesTrace.add(ranks[i].getPosition());
			} else {
				tiesTrace.add(ranks[i].getPosition());
			}

			out[ranks[i].getPosition()] = (double) pos;
		}

		if (tiesTrace.size() > 1) {
			this.resolveTie(out, tiesTrace);
		}

		if (this.nanStrategy == NaNStrategy.FIXED) {
			this.restoreNaNs(out, nanPositions);
		}

		return out;
	}

	private IntDoublePair[] removeNaNs(IntDoublePair[] ranks) {
		if (!this.containsNaNs(ranks)) {
			return ranks;
		} else {
			IntDoublePair[] outRanks = new IntDoublePair[ranks.length];
			int j = 0;

			for (int i = 0; i < ranks.length; ++i) {
				if (Double.isNaN(ranks[i].getValue())) {
					for (int k = i + 1; k < ranks.length; ++k) {
						ranks[k] = new IntDoublePair(ranks[k].getValue(), ranks[k].getPosition() - 1);
					}
				} else {
					outRanks[j] = new IntDoublePair(ranks[i].getValue(), ranks[i].getPosition());
					++j;
				}
			}

			IntDoublePair[] returnRanks = new IntDoublePair[j];
			System.arraycopy(outRanks, 0, returnRanks, 0, j);
			return returnRanks;
		}
	}

	private void recodeNaNs(IntDoublePair[] ranks, double value) {
		IntStream.range(0, ranks.length).parallel().forEach(i -> {
			if (Double.isNaN(ranks[i].getValue())) {
				ranks[i] = new IntDoublePair(value, ranks[i].getPosition());
			}
		});
	}

	private boolean containsNaNs(IntDoublePair[] ranks) {
		for (int i = 0; i < ranks.length; ++i) {
			if (Double.isNaN(ranks[i].getValue())) {
				return true;
			}
		}
		return false;
	}

	private void resolveTie(double[] ranks, List<Integer> tiesTrace) {
		double c = ranks[(Integer) tiesTrace.get(0)];
		int length = tiesTrace.size();
		Iterator iterator;
		long f;
		switch (this.tiesStrategy) {
			case AVERAGE:
				this.fill(ranks, tiesTrace, (2.0D * c + (double) length - 1.0D) / 2.0D);
				break;
			case MAXIMUM:
				this.fill(ranks, tiesTrace, c + (double) length - 1.0D);
				break;
			case MINIMUM:
				this.fill(ranks, tiesTrace, c);
				break;
			case RANDOM:
				iterator = tiesTrace.iterator();

				for (f = FastMath.round(c); iterator.hasNext(); ranks[(Integer) iterator.next()] = (double) this.randomData.nextLong(f, f + (long) length - 1L)) {
				}

				return;
			case SEQUENTIAL:
				iterator = tiesTrace.iterator();
				f = FastMath.round(c);

				for (int var9 = 0; iterator.hasNext(); ranks[(Integer) iterator.next()] = (double) (f + (long) (var9++))) {
				}

				return;
			default:
				throw new MathInternalError();
		}

	}

	private void fill(double[] data, List<Integer> tiesTrace, double value) {
		for (Iterator iterator = tiesTrace.iterator(); iterator.hasNext(); data[(Integer) iterator.next()] = value) {
		}

	}

	private void restoreNaNs(double[] ranks, List<Integer> nanPositions) {
		if (nanPositions.size() != 0) {
			for (Iterator iterator = nanPositions.iterator(); iterator.hasNext(); ranks[(Integer) iterator.next()] = 0.0D / 0.0) {
			}

		}
	}

	private List<Integer> getNanPositions(IntDoublePair[] ranks) {
		ArrayList<Integer> out = new ArrayList();

		for (int i = 0; i < ranks.length; ++i) {
			if (Double.isNaN(ranks[i].getValue())) {
				out.add(i);
			}
		}

		return out;
	}

	static {
		DEFAULT_NAN_STRATEGY = NaNStrategy.FAILED;
		DEFAULT_TIES_STRATEGY = TiesStrategy.AVERAGE;
	}

	private static class IntDoublePair implements Comparable<IntDoublePair> {
		private final double value;
		private final int position;

		IntDoublePair(double value, int position) {
			this.value = value;
			this.position = position;
		}

		public int compareTo(IntDoublePair other) {
			return Double.compare(this.value, other.value);
		}

		public double getValue() {
			return this.value;
		}

		public int getPosition() {
			return this.position;
		}
	}


}
