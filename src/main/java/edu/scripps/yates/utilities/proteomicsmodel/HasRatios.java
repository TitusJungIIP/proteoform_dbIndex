package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.Set;

public interface HasRatios {
	public Set<Ratio> getRatios();

	public boolean addRatio(Ratio ratio);
}
