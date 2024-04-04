package edu.scripps.yates.utilities.proteomicsmodel.factories;

import java.io.Serializable;

import org.apache.commons.lang.builder.HashCodeBuilder;

import edu.scripps.yates.utilities.proteomicsmodel.Condition;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.enums.CombinationType;

/**
 * This class represents a ratio of the expression of a {@link Protein} between
 * two {@link ConditionEx}s
 *
 * @author Salva
 *
 */
public class RatioEx implements Ratio, Serializable {

	/**
	 *
	 */
	private static final long serialVersionUID = -4788750494296885167L;
	private final double value;
	private Condition condition1;
	private Condition condition2;
	private final String description;
	private Score score;
	private CombinationType combinationType;
	private final AggregationLevel aggregationLevel;
	private double stdev;

	public RatioEx(double value, Condition condition1, Condition condition2, String ratioDescription,
			AggregationLevel aggregationLevel) {
		super();
		this.value = value;
		this.condition1 = condition1;
		this.condition2 = condition2;
		description = ratioDescription;
		this.aggregationLevel = aggregationLevel;
	}

	public RatioEx(double value, Condition condition1, Condition condition2, CombinationType combinationType,
			String ratioDescription, AggregationLevel aggregationLevel) {
		super();
		this.value = value;
		this.condition1 = condition1;
		this.condition2 = condition2;
		description = ratioDescription;
		this.combinationType = combinationType;
		this.aggregationLevel = aggregationLevel;
	}

	/**
	 * @return the value
	 */
	@Override
	public double getValue() {
		return value;
	}

	/**
	 * @return the condition1
	 */
	@Override
	public Condition getCondition1() {
		return condition1;
	}

	/**
	 * @return the condition2
	 */
	@Override
	public Condition getCondition2() {
		return condition2;
	}

	/**
	 * @return the description
	 */
	@Override
	public String getDescription() {
		return description;
	}

	@Override
	public Score getAssociatedConfidenceScore() {
		return score;
	}

	public void setAssociatedConfidenceScore(Score score) {
		this.score = score;
	}

	public void setStandardDeviationOfLog2Ratios(double stdev) {
		this.stdev = stdev;
	}

	public void setStdev(double stdev) {
		this.stdev = stdev;
	}

	@Override
	public CombinationType getCombinationType() {
		return combinationType;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object arg0) {
		if (arg0 instanceof Ratio) {
			final Ratio ratio = (Ratio) arg0;
			if (ratio.getValue() != getValue())
				return false;
			if (!ratio.getCondition1().equals(getCondition1()))
				return false;
			if (!ratio.getCondition2().equals(getCondition2()))
				return false;
			if (ratio.getDescription() != null && getDescription() != null
					&& !ratio.getDescription().equals(getDescription()))
				return false;
			return true;
		}

		return super.equals(arg0);

	}

	@Override
	public int hashCode() {
		final HashCodeBuilder hashCodeBuilder = new HashCodeBuilder();

		hashCodeBuilder.append(getValue());
		if (getCondition1() != null)
			hashCodeBuilder.append(getCondition1().hashCode());
		if (getCondition2() != null)
			hashCodeBuilder.append(getCondition2().hashCode());
		if (getDescription() != null) {
			hashCodeBuilder.append(getDescription());
		}

		return hashCodeBuilder.toHashCode();

	}

	@Override
	public AggregationLevel getAggregationLevel() {
		return aggregationLevel;
	}

	/**
	 * @return the score
	 */
	public Score getScore() {
		return score;
	}

	/**
	 * @param score the score to set
	 */
	public void setScore(Score score) {
		this.score = score;
	}

	/**
	 * @param combinationType the combinationType to set
	 */
	public void setCombinationType(CombinationType combinationType) {
		this.combinationType = combinationType;
	}

	/**
	 * @param condition1 the condition1 to set
	 */
	public void setCondition1(Condition condition1) {
		this.condition1 = condition1;
	}

	/**
	 * @param condition2 the condition2 to set
	 */
	public void setCondition2(Condition condition2) {
		this.condition2 = condition2;
	}

	@Override
	public double getStandardDeviationOfLog2Ratios() {
		return stdev;
	}

}
