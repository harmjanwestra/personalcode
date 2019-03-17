package nl.harmjanwestra.playground.legacy.vcf.filter.genotypefilters;

import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;

/**
 * Created by hwestra on 2/8/16.
 */
public interface VCFGenotypeFilter {
	void filter(VCFVariant variant);
}
