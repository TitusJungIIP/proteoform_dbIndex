package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.Set;

public interface HasMsRuns {
	public Set<MSRun> getMSRuns();

	public boolean addMSRun(MSRun msRun);
}
