package nl.harmjanwestra.playground.conv;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

public class ConvertVCFToTTTest {

	public static void main(String[] args) {
		String testfile = "D:\\tmp\\chr22.dose.vcf.gz";
		String testout = "D:\\tmp\\uot\\";

		ConvertVCFToTT v = new ConvertVCFToTT();
		String[] allowedVQSROrFilter = new String[]{
				"PASS", "GENOTYPED",
				"GENOTYPED_ONLY",
				"VQSRTrancheINDEL99.90to99.95",
				"VQSRTrancheINDEL99.00to100.00+",
				"VQSRTrancheINDEL99.00to100.00",
				"VQSRTrancheINDEL99.95to100.00+",
				"VQSRTrancheINDEL99.95to100.00",
				"VQSRTrancheSNP99.80to100.00",
				"VQSRTrancheSNP99.80to99.90",
				"VQSRTrancheSNP99.90to99.95",
				"VQSRTrancheSNP99.95to100.00+",
				"VQSRTrancheSNP99.95to100.00",
				"VQSRTrancheSNP99.80to100.00",
				"VQSRTrancheSNP99.80to100.00+",
				"VQSRTrancheSNP99.80to100.00"
		};
		HashSet<String> allowedFilter = new HashSet<String>();
		allowedFilter.addAll(Arrays.asList(allowedVQSROrFilter));
		try {
			v.run(testfile, testout, true, true, 0.01, 0.9, 0, 0, 1, 20, allowedFilter, -0.3, 0.3);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
