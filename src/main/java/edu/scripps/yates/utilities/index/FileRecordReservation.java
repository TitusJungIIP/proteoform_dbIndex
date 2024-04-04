package edu.scripps.yates.utilities.index;

import java.io.File;

public class FileRecordReservation {

	private long currentLastPosition;

	public FileRecordReservation(File file) {
		this.currentLastPosition = file.length();
	}

	/**
	 * It returns the position in which the record should be written in the file
	 * 
	 * @param recordToAppend
	 * @return
	 */
	public synchronized long reserveRecord(byte[] recordToAppend) {
		long ret = currentLastPosition;
		currentLastPosition += recordToAppend.length;
		return ret;
	}

	public synchronized long getCurrentposition() {
		return currentLastPosition;
	}
}
