package edu.scripps.yates.utilities.progresscounter;

import java.text.DecimalFormat;

import org.apache.commons.lang.mutable.MutableLong;
import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.dates.DatesUtil;

public class ProgressCounter {
	private final static Logger log = Logger.getLogger(ProgressCounter.class);
	private long total;
	private MutableLong count;
	private String previousPercentage = "";
	private final ProgressPrintingType progressPrintingType;
	private final DecimalFormat df;
	private boolean showRemainingTime;
	private long t1 = -1;
	private ProgressNumberFormatter formatter;
	private String suffix;
	private String percentageToAssign;

	/**
	 * Progress counter showing time remaining by default
	 * 
	 * @param total
	 * @param progressPrintingType
	 * @param numDecimals
	 */
	public ProgressCounter(long total, ProgressPrintingType progressPrintingType, int numDecimals) {
		this(total, progressPrintingType, numDecimals, true);
	}

	/**
	 * Progress counter with all parameters to set.
	 * 
	 * @param total
	 * @param progressPrintingType
	 * @param numDecimals
	 * @param showTimeRemaining
	 */
	public ProgressCounter(long total, ProgressPrintingType progressPrintingType, int numDecimals,
			boolean showTimeRemaining) {
		this.total = total;
		count = new MutableLong(0);
		this.progressPrintingType = progressPrintingType;
		final StringBuilder sb = new StringBuilder("#");
		for (int i = 0; i < numDecimals; i++) {
			if ("#".equals(sb.toString())) {
				sb.append(".");
			}
			sb.append("#");
		}
		df = new DecimalFormat(sb.toString());
		showRemainingTime = showTimeRemaining;
	}

	public boolean isShowRemainingTime() {
		return showRemainingTime;
	}

	public void setShowRemainingTime(boolean showRemainingTime) {
		this.showRemainingTime = showRemainingTime;
	}

	/**
	 * Increments the counter by increment
	 * 
	 * @param increment
	 */
	public void increment(long increment) {
		previousPercentage = percentageToAssign;
		count.add(increment);
		startTimeWithFirstIncrement();

	}

	public void setProgress(long progress) {
		previousPercentage = percentageToAssign;
		count.setValue(progress);
		startTimeWithFirstIncrement();

	}

	/**
	 * increments the counter by 1
	 */
	public void increment() {
		previousPercentage = percentageToAssign;
		count.add(1);
		startTimeWithFirstIncrement();

	}

	private void startTimeWithFirstIncrement() {
		// start time with the first increment
		if (showRemainingTime && t1 == -1) {
			log.debug("taking time at the first increment");
			t1 = System.currentTimeMillis();
		}
	}

	private String getRemainingTime() {
		if (showRemainingTime) {
			final long currentTimeMillis = System.currentTimeMillis();
			final long timeConsumed = currentTimeMillis - t1;
			final double timeConsumedPerItem = 1.0 * timeConsumed / getCount();
			final double estimatedRemainingTime = timeConsumedPerItem * (1.0 * total - getCount());
			return DatesUtil.getDescriptiveTimeFromMillisecs(estimatedRemainingTime);
		}
		return "";
	}

	/**
	 * Returns the progress string or empty string according to the
	 * ProgressPrintingType. If the {@link ProgressPrintingType} is
	 * PERCENTAGE_STEPS, then this function will only return the progress string
	 * every time the percentage changes. The percentage changes according to the
	 * number of decimals stated in the constructor.
	 * 
	 * @return
	 */
	public String printIfNecessary() {

		final StringBuilder sb = new StringBuilder();
		if (progressPrintingType == ProgressPrintingType.EVERY_STEP) {
			final double percentage = getPercentage();
			sb.append(getFormatter().format(count.longValue())).append("/").append(getFormatter().format(total))
					.append(" (").append(df.format(percentage)).append("%)").toString();
			if (showRemainingTime) {
				sb.append(" (" + getRemainingTime() + " remaining...)");
			}
			if (suffix != null) {
				sb.append(" ").append(suffix);
			}
			return sb.toString();
		} else if (progressPrintingType == ProgressPrintingType.PERCENTAGE_STEPS) {
			final String percentage = df.format(getPercentage());
			if (!percentage.equals(previousPercentage)) {
				sb.append(getFormatter().format(count.longValue())).append("/").append(getFormatter().format(total))
						.append(" (").append(percentage).append("%)").toString();
				if (showRemainingTime) {
					sb.append(" (" + getRemainingTime() + " remaining...)");
				}
				percentageToAssign = percentage;
				if (suffix != null) {
					sb.append(" ").append(suffix);
				}
				return sb.toString();
			}
		}
		return "";
	}

	public long getCount() {
		return count.longValue();
	}

	public double getPercentage() {
		final double d = count.longValue() * 100.0 / total;
		return d;
	}

	public void setTotal(long max) {
		total = max;
		count = new MutableLong(0);
		t1 = -1;
	}

	public void setCount(long count2) {
		previousPercentage = percentageToAssign;
		count.setValue(count2);
		if (count2 > 0) {
			startTimeWithFirstIncrement();
		}
	}

	public void addCount(long count2) {
		previousPercentage = percentageToAssign;
		count.add(count2);
		startTimeWithFirstIncrement();
	}

	public void setProgressNumberFormatter(ProgressNumberFormatter formatter) {
		this.formatter = formatter;
	}

	public ProgressNumberFormatter getFormatter() {
		if (formatter == null) {
			formatter = new ProgressNumberFormatter() {

				@Override
				public String format(long number) {
					return String.valueOf(number);
				}
			};
		}
		return formatter;
	}

	public String getSuffix() {
		return suffix;
	}

	public void setSuffix(String suffix) {
		this.suffix = suffix;
	}
}
