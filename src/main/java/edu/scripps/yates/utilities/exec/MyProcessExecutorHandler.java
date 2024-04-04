package edu.scripps.yates.utilities.exec;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

public class MyProcessExecutorHandler implements ProcessExecutorHandler {
	private final static Logger log = Logger.getLogger(MyProcessExecutorHandler.class);

	private final List<String> standardOutputMessages = new ArrayList<String>();
	private final List<String> standardErrorMessages = new ArrayList<String>();
	private final static String[] knownErrorMessagesStarts = { "error:", "indexerror:", "list index out of range" };

	@Override
	public void onStandardOutput(String msg) {
		log.info(msg);
		standardOutputMessages.add(msg);
	}

	@Override
	public void onStandardError(String msg) {
		log.error(msg);
		standardErrorMessages.add(msg);
	}

	public List<String> getStandardOutputMessages() {
		return standardOutputMessages;
	}

	public List<String> getStandardErrorMessages() {
		return standardErrorMessages;
	}

	public boolean containsErrorMessage() {
		for (final String msg : standardErrorMessages) {
			for (final String errorMsg : knownErrorMessagesStarts) {
				if (msg.toLowerCase().contains(errorMsg)) {
					return true;
				}
			}

		}
		for (final String msg : standardOutputMessages) {
			for (final String errorMsg : knownErrorMessagesStarts) {
				if (msg.toLowerCase().contains(errorMsg)) {
					return true;
				}
			}
		}
		return false;
	}

	public String getErrorMessage() {
		for (final String msg : standardErrorMessages) {
			if (msg.toLowerCase().contains("error:")) {
				return msg;
			}
		}
		for (final String msg : standardOutputMessages) {
			if (msg.toLowerCase().contains("error:")) {
				return msg;
			}
		}
		return null;
	}

	public void resetMesssages() {
		standardErrorMessages.clear();
		standardOutputMessages.clear();
	}
}
