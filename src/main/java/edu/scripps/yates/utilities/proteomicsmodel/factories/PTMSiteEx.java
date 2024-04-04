package edu.scripps.yates.utilities.proteomicsmodel.factories;

import java.io.Serializable;

import edu.scripps.yates.utilities.proteomicsmodel.PTMPosition;
import edu.scripps.yates.utilities.proteomicsmodel.PTMSite;
import edu.scripps.yates.utilities.proteomicsmodel.Score;

public class PTMSiteEx implements PTMSite, Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -4135959343681112419L;
	private final String aa;
	private final int position;
	private Score score;
	private PTMPosition ptmPosition;

	public PTMSiteEx(String aa, int position, PTMPosition ptmPosition) {
		this.aa = aa;
		this.position = position;
		this.ptmPosition = ptmPosition;
	}

	@Override
	public String getAA() {
		return aa;
	}

	@Override
	public int getPosition() {
		return position;
	}

	@Override
	public Score getScore() {
		return score;
	}

	/**
	 * @param score the score to set
	 */
	public void setScore(Score score) {
		this.score = score;
	}

	@Override
	public PTMPosition getPTMPosition() {
		return ptmPosition;
	}

	@Override
	public void setPTMPosition(PTMPosition ptmPosition) {
		this.ptmPosition = ptmPosition;
	}

}
