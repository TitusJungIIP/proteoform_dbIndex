package edu.scripps.yates.utilities.exec;

import java.io.IOException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import org.apache.commons.exec.CommandLine;
import org.apache.log4j.Logger;

public class CommandLineRunner {
	private final static Logger log = Logger.getLogger(CommandLine.class);
	private final MyProcessExecutorHandler handler = new MyProcessExecutorHandler();
	private Long processExitCode;

	public void runCommand(CommandLine commandLine, long timeout)
			throws IOException, InterruptedException, ExecutionException {
		final String commandString = commandLine.toString();
		log.info("Running: " + commandString);
		handler.resetMesssages();
		final Future<Long> runProcess = ProcessExecutor.runProcess(commandLine, handler, timeout);
		while (!runProcess.isDone() && !runProcess.isCancelled()) {
			Thread.sleep(1000);
		}
		processExitCode = runProcess.get();
		if (handler.containsErrorMessage()) {
			processExitCode = -1l;
		}
		log.info("Process exitValue: " + processExitCode);
	}

	public String getErrorMessage() {
		return handler.getErrorMessage();
	}

	public boolean containsErrorMessage() {
		return handler.containsErrorMessage();
	}

	public Long getProcessExitCode() {
		return processExitCode;
	}
}
