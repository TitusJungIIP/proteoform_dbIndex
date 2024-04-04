package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.Set;

public interface HasScores {
	public Set<Score> getScores();

	public boolean addScore(Score score);

}
