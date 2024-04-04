package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.List;

public interface HasPTMs {
	public List<PTM> getPTMs();

	public boolean addPTM(PTM newPtm);
}
