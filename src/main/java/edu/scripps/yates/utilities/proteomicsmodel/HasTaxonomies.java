package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.Set;

public interface HasTaxonomies {
	public boolean addTaxonomy(String taxonomy);

	public Set<String> getTaxonomies();

	public boolean isIgnoreTaxonomy();

	public void setIgnoreTaxonomy(boolean ignoreTaxonomy);
}
