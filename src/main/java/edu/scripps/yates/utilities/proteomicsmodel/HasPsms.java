package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.List;

public interface HasPsms {
	public List<PSM> getPSMs();

	public boolean addPSM(PSM psm, boolean recursively);

	public Integer getSpectrumCount();

	public void setSpectrumCount(Integer spc);
}
