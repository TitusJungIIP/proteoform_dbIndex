package edu.scripps.yates.utilities.memory;

import java.text.DecimalFormat;

import edu.scripps.yates.utilities.files.FileUtils;

public class MemoryUsageReport {
	private final static DecimalFormat df = new DecimalFormat("#.0");
	public static final double RECOMMENDED_MINIMUM_MEMORY_PERCENTAGE = 10.0;

	public static String getMemoryUsageReport() {

		final Runtime runtime = Runtime.getRuntime();
		final long totalMemory = runtime.totalMemory();
		final long maxMemory = runtime.maxMemory();
		final long freeMemory = runtime.freeMemory();
		// long freeMemory = runtime.freeMemory();

		final long usedMemory = totalMemory - freeMemory;
		final long realFreeMemory = maxMemory - usedMemory;
		final double usedPercentage = usedMemory * 100.0 / maxMemory;
		final double realFreePercentage = realFreeMemory * 100.0 / maxMemory;
		return FileUtils.getDescriptiveSizeFromBytes(usedMemory) + " used (" + df.format(usedPercentage) + "%), "
				+ FileUtils.getDescriptiveSizeFromBytes(realFreeMemory) + " free (" + df.format(realFreePercentage)
				+ "%)";
	}

	public static String getUsedMemoryDescriptiveString() {
		final Runtime runtime = Runtime.getRuntime();
		final long totalMemory = runtime.totalMemory();
		final long maxMemory = runtime.maxMemory();
		final long freeMemory = runtime.freeMemory();
		final long usedMemory = totalMemory - freeMemory;
		final long realFreeMemory = maxMemory - usedMemory;
		return FileUtils.getDescriptiveSizeFromBytes(realFreeMemory);
	}

	public static String getFreeMemoryDescriptiveString() {
		final Runtime runtime = Runtime.getRuntime();
		final long totalMemory = runtime.totalMemory();
		final long freeMemory = runtime.freeMemory();
		final long usedMemory = totalMemory - freeMemory;
		return FileUtils.getDescriptiveSizeFromBytes(usedMemory);
	}

	public static double getUsedMemoryPercentage() {
		final Runtime runtime = Runtime.getRuntime();
		final long totalMemory = runtime.totalMemory();
		final long maxMemory = runtime.maxMemory();
		final long freeMemory = runtime.freeMemory();

		final long usedMemory = totalMemory - freeMemory;
		final double usedPercentage = usedMemory * 100.0 / maxMemory;
		return usedPercentage;
	}

	public static double getFreeMemoryPercentage() {
		final Runtime runtime = Runtime.getRuntime();
		final long totalMemory = runtime.totalMemory();
		final long maxMemory = runtime.maxMemory();
		final long freeMemory = runtime.freeMemory();

		final long usedMemory = totalMemory - freeMemory;
		final long realFreeMemory = maxMemory - usedMemory;
		final double freePercentage = realFreeMemory * 100.0 / maxMemory;
		return freePercentage;
	}
}
