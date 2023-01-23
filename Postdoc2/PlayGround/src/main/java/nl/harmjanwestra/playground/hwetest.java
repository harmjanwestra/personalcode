package nl.harmjanwestra.playground;

import htsjdk.tribble.util.popgen.HardyWeinbergCalculation;
import umcg.genetica.math.stats.HWE;

public class hwetest {

	public static void main(String[] args) {
		int aa = 1808;
		int ab = 86;
		int bb = 0;

		// nrhets: 86, nrhomA: 1808, nrhomB: 0

		double p = HardyWeinbergCalculation.hwCalculate(aa, ab, bb);
		double p2 = HWE.calculateExactHWEPValue(ab, aa, bb);
		double p3 = hwetest.calculateExactHWEPValue(ab, aa, bb);


		/*
		1       1:14590:.:G,A   G       A       1808    86      0       0.0454065       0.0443757       0.623593
1:14590:.       FailedVariantStats      CR=1.0;MAF=0.023;HWE=0.000e+00
		 */

		System.out.println(p);
		System.out.println(p2);
		System.out.println(p3);
	}


	public static double calculateExactHWEPValue(int obs_hets, int obs_hom1, int obs_hom2) {
		//System.out.println("Starting exact HWE:\t" + obs_hets + "\t" + obs_hom1 + "\t" + obs_hom2);

		int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
		int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

		int rare_copies = 2 * obs_homr + obs_hets;
		int l_genotypes = obs_hets + obs_homc + obs_homr;

		if (l_genotypes == 0) {
			return -1;
		}

		System.out.println(rare_copies + "\t" + l_genotypes);

		double[] het_probs = new double[rare_copies + 1];
		int i;

		/* start at midpoint */
		int mid = rare_copies * (2 * l_genotypes - rare_copies) / (2 * l_genotypes);
		/* check to ensure that midpoint and rare alleles have same parity */
		if (mid % 2 != rare_copies % 2) {
			mid++;
		}

		System.out.println(mid);
		int curr_hets = mid;
		int curr_homr = (rare_copies - mid) / 2;
		int curr_homc = l_genotypes - curr_hets - curr_homr;

		System.out.println(curr_hets + "\t" + curr_homr + "\t" + curr_homc);

		het_probs[mid] = 1.0;
		double sum = het_probs[mid];
		for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
			het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
			System.out.println(curr_hets + "\t" + het_probs[curr_hets - 2]);
			sum += het_probs[curr_hets - 2];
			/* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
			curr_homr++;
			curr_homc++;
		}

		curr_hets = mid;
		curr_homr = (rare_copies - mid) / 2;
		curr_homc = l_genotypes - curr_hets - curr_homr;
		System.out.println(curr_hets + "\t" + curr_homr + "\t" + curr_homc);
		for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
			het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0));
			System.out.println(curr_hets + "\t" + het_probs[curr_hets - 2]);
			sum += het_probs[curr_hets + 2];
			curr_homr--;
			curr_homc--;
		}

		for (i = 0; i <= rare_copies; i++) {
			het_probs[i] /= sum;
			System.out.println(i + "\t" + het_probs[i] + "\t" + sum);
		}

		double p_hwe = 0.0;
		for (i = 0; i <= rare_copies; i++) {
			if (het_probs[i] <= het_probs[obs_hets]) {
				p_hwe += het_probs[i];
			}
		}

		p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

		return p_hwe;
	}
}
