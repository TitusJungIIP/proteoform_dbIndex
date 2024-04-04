package edu.scripps.yates.utilities.swing;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * This abstract class defines the implementor as a class that can be used by
 * {@link CommandLineProgramGuiEnclosable} to generate a GUI from a command line
 * program
 * 
 * @author salvador
 *
 */
public abstract class CommandLineProgramGuiEnclosable<T> implements StatusListener<T> {
	protected static final String GUI = "gui";
	protected static final String GUI_LONG = "graphical_interface";

	private Options options;
	private boolean readyForRun;
	protected AutomaticGUICreator gui;
	private CommandLine commandLine;
	private final boolean usingGUI;

	/**
	 * 
	 * @param mainArgs the arguments passed to the
	 *                 {@code public static void main(String[] args) } method
	 * @throws ParseException
	 * @throws SomeErrorInParametersOcurred
	 */
	public CommandLineProgramGuiEnclosable(String[] mainArgs) throws ParseException, SomeErrorInParametersOcurred {
		usingGUI = hasGUIOption(mainArgs);
		try {
			final CommandLineParser parser = new DefaultParser();
			final List<Option> optionList = defineCommandLineOptions();
			options = new Options();
			// add new option gui
			options.addOption(GUI, GUI_LONG, false,
					"Launch this tool with a graphical interface. All other command line parameters will be ignored.");

			if (usingGUI) {

				for (final Option option : optionList) {
					final Option optionalOption = (Option) option.clone();
					optionalOption.setRequired(false);
					options.addOption(optionalOption);
				}
				mainArgs = removeGUIOptionFromArguments(mainArgs);
				parser.parse(options, mainArgs);
				options = new Options();
				for (final Option option : optionList) {
					options.addOption(option);
				}
				startGUI();
				// set to false
				readyForRun = false;
			} else {
				readyForRun = true;
				for (final Option option : optionList) {
					if (!option.hasArg()) {
						option.setRequired(false);
					}
					options.addOption(option);
				}

				mainArgs = removeGUIOptionFromArguments(mainArgs);
				commandLine = parser.parse(options, mainArgs);
				initTool(commandLine);
			}
		} catch (final MissingOptionException e) {
			if (!usingGUI) {
				errorInParameters(e.getMessage() + ".  ***You can also run this tool with '-" + GUI
						+ "' option to launch the graphical interface.***");
			} else {
				throw e;
			}
		}
	}

	public CommandLine getCommandLine() {
		return commandLine;
	}

	private boolean hasGUIOption(String[] mainArgs) {
		for (final String arg : mainArgs) {
			if (arg.startsWith("-")) {
				if (arg.equals("-" + GUI)) {
					return true;
				}
			} else if (arg.startsWith("--")) {
				if (arg.equals("--" + GUI_LONG)) {
					return true;
				}
			}
		}
		return false;

	}

	private String[] removeGUIOptionFromArguments(String[] mainArgs) {
		final List<String> ret = new ArrayList<String>();
		for (final String arg : mainArgs) {
			if (arg.startsWith("-")) {
				if (!arg.equals("-" + GUI)) {
					ret.add(arg);
				}
			} else if (arg.startsWith("--")) {
				if (!arg.equals("--" + GUI_LONG)) {
					ret.add(arg);
				}
			} else {
				ret.add(arg);
			}
		}
		return ret.toArray(new String[0]);
	}

	private void startGUI() {
		gui = new AutomaticGUICreator(this);
		gui.setVisible(true);

	}

	@Override
	public void onStatusUpdate(String statusMessage) {
		if (gui != null) {
			gui.getStatusPrintStream().println(statusMessage);
		}

	}

	@Override
	public void onStatusUpdate(String statusMessage, Object param) {
		if (gui != null) {
			gui.getStatusPrintStream().println(statusMessage);
		}

	}

	protected abstract List<Option> defineCommandLineOptions();

	protected Options getCommandLineOptions() {
		return options;
	};

	/**
	 * Method invoked by the graphical interface when the user press the RUN button.
	 * Do NOT invoke this method programmatically. Use <code>safeRun()</code>
	 * instead.
	 * 
	 * @throws Exception
	 */
	public abstract void run();

	/**
	 * Title of the frame window
	 * 
	 * @return
	 */
	public abstract String getTitleForFrame();

	/**
	 * Init the tool with the options set in the command line.<br>
	 * Use this method to further check the validity of the argument captured from
	 * the GUI (such as checking whether a file exists or not), and it is not valid,
	 * call to
	 * 
	 * <pre>
	 * {@code errorInParameters(error_message);}
	 * </pre>
	 * 
	 * For example:
	 * 
	 * <pre>
	 * {@code 
	* if (cmd.hasOption("in")) { 
	*    final String inputFilePath = cmd.getOptionValue("in"); 
	* } else { 
	* 	 super.errorInParameters("Input file path is is missing"); 
	* }
	 * </pre>
	 * 
	 * Do not use any final variable in this method since this will be executed
	 * before the initialization of the class (it is called in the super
	 * constructor). @param cmd command line object
	 * 
	 * @throws SomeErrorInParametersOcurred
	 */
	protected abstract void initToolFromCommandLineOptions(CommandLine cmd) throws SomeErrorInParametersOcurred;

	void initTool(CommandLine cmd) throws SomeErrorInParametersOcurred {
		commandLine = cmd;
		initToolFromCommandLineOptions(cmd);
	}

	/**
	 * This method should print the usage of the command line such as:<br>
	 * programName -i [input file] -an [value] ..."
	 * 
	 * @return
	 */
	public abstract String printCommandLineSintax();

	/**
	 * Override this method if you want something different than:<br>
	 * "Contact Salvador Martinez-Bartolome at salvador@scripps.edu for more help"
	 * 
	 * @return
	 */
	public String printFooter() {
		return "Contact Salvador Martinez-Bartolome at salvador@scripps.edu for more help";
	}

	/**
	 * Prints an error regarding the parameter and subsequently throws
	 * {@link SomeErrorInParametersOcurred}
	 * 
	 * @param header
	 * @throws SomeErrorInParametersOcurred
	 */
	protected void errorInParameters(String header) throws SomeErrorInParametersOcurred {
		// automatically generate the help statement
		final HelpFormatter formatter = new HelpFormatter();
		String printCommandLineSintax = printCommandLineSintax();
		final boolean printAutomaticUsage = printCommandLineSintax == null;
		if (printCommandLineSintax == null) {
			printCommandLineSintax = "java -jar *.jar";
		}
//		if (header == null) {
//			formatter.printHelp(250,printCommandLineSintax, options);
//		} else {
//		if (header != null) {
//			header = "\n\n************\n" + header + "\n************\n\n";
//		}
		formatter.printHelp(250, printCommandLineSintax, header, options, printFooter(), printAutomaticUsage);
//		}
		throw new SomeErrorInParametersOcurred(header);
	}

	/**
	 * This method should be the one to invoke to start the process from an invoker
	 * other than the GUI rather than run()
	 * 
	 * @throws Exception
	 */
	public void safeRun() throws DoNotInvokeRunMethod {
		if (!readyForRun) {
			throw new DoNotInvokeRunMethod();
		} else {
			run();
			System.exit(0);
		}

	}

	public boolean isUsingGUI() {
		return usingGUI;
	}
}
