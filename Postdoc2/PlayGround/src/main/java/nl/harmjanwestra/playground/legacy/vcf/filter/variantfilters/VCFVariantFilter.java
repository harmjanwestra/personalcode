package nl.harmjanwestra.playground.legacy.vcf.filter.variantfilters;

import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;

/**
 * Created by hwestra on 3/21/17.
 */
public interface VCFVariantFilter {
	boolean passesThreshold(VCFVariant variant);
}
