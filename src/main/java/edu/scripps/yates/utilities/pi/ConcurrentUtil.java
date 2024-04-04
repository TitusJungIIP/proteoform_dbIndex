package edu.scripps.yates.utilities.pi;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.TimeUnit;

/**
 * Class for managing multi-threading
 * 
 * @author salvador
 *
 */
public class ConcurrentUtil {
	/**
	 * This method will call shutdown and will wait for the threads of the
	 * {@link ExecutorService} to finish for the time stated in the
	 * awaitTerminationMilliseconds parameter. If they are not finished by that
	 * time, the method shutdownNow will be called to force all threads to be
	 * interrupted.
	 * 
	 * @param executor
	 * @param awaitTerminationMilliseconds
	 */
	public static void stop(ExecutorService executor, long awaitTerminationMilliseconds) {
		try {
			executor.shutdown();
			executor.awaitTermination(awaitTerminationMilliseconds, TimeUnit.MILLISECONDS);
		} catch (final InterruptedException e) {
			System.err.println("termination interrupted");
		} finally {
			if (!executor.isTerminated()) {
				System.err.println("killing non-finished tasks");
			}
			executor.shutdownNow();
		}
	}

	/**
	 * This will call to sleep to the current thread and the
	 * {@link InterruptedException} will be catch and a {@link RuntimeException}
	 * will wrap it and will be thrown
	 * 
	 * @param miliseconds
	 */
	public static void sleep(long miliseconds) {
		try {
			Thread.sleep(miliseconds);
		} catch (final InterruptedException e) {
			throw new RuntimeException(e);
		}
	}

}
