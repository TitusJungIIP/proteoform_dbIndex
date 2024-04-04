package edu.scripps.yates.utilities.swing;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Label;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.List;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.WindowConstants;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.Logger;
import org.springframework.core.io.ClassPathResource;

import edu.scripps.yates.utilities.appversion.AppVersion;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.properties.PropertiesUtil;
import gnu.trove.map.TMap;
import gnu.trove.map.hash.THashMap;

public class AutomaticGUICreator extends JFrame {

	/**
	 * 
	 */
	private static final long serialVersionUID = 4829866413822345381L;
	private static Logger log = Logger.getLogger(AutomaticGUICreator.class);
	private static AppVersion version;
	private final CommandLineProgramGuiEnclosable program;
	private final TMap<String, JComponent> componentsByOption = new THashMap<String, JComponent>();
	private final JTextArea status;
	private final JSplitPane splitPane;
	private final MyPrintStream statusPrintStream;
	private static final SimpleDateFormat format = new SimpleDateFormat("yyy-MM-dd HH:mm:ss,SSS");

	public AutomaticGUICreator(CommandLineProgramGuiEnclosable program) {
		super();
		final AppVersion version2 = getVersion();
		String title = program.getTitleForFrame();
		if (version2 != null) {
			title += " - v" + version2;
		}
		super.setTitle(title);
		setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		this.program = program;

		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (UnsupportedLookAndFeelException | ClassNotFoundException | InstantiationException
				| IllegalAccessException e) {
			e.printStackTrace();
		}
		getContentPane().setLayout(new BorderLayout(10, 10));
		final JPanel headerPanel = new JPanel();
		headerPanel.add(new Label("Graphical Interface for " + program.getTitleForFrame()));
		getContentPane().add(headerPanel, BorderLayout.NORTH);

		final GridBagLayout layout = new GridBagLayout();
		layout.columnWeights = new double[] { 0.1, 0.1, 0.5, 0.3 };
		layout.columnWidths = new int[] { 30, 20, 250, 100 };
		final JPanel componentsPanel = new JPanel(layout);
		componentsPanel.setBorder(BorderFactory.createTitledBorder("Parameters:"));
		final JScrollPane scroll2 = new JScrollPane(new MyScrollablePanel(componentsPanel));
		final Options options = program.getCommandLineOptions();
		final Collection<Option> optionList = options.getOptions();
		int y = 0;
		for (final Option option : optionList) {
			if (option.getOpt().equals(CommandLineProgramGuiEnclosable.GUI)) {
				// ignore it
				continue;
			}
			final boolean required = option.isRequired();
			final JLabel label = getLabelForOptionName(option);
			componentsPanel.add(label, getGridBagConstraints(0, y, GridBagConstraints.NORTHEAST));

			final JLabel labelRequired = new JLabel();
			if (required) {
				labelRequired.setText("[Required]");
				labelRequired.setToolTipText("This parameter is required");
			} else {
				labelRequired.setText("[Optional]");
				labelRequired.setToolTipText("This parameter is optional. The tool can run with this parameter empty");
			}
			componentsPanel.add(labelRequired, getGridBagConstraints(1, y, GridBagConstraints.NORTHWEST));

			final JComponent component = getComponentForOption(option, componentsByOption);
			componentsPanel.add(component, getGridBagConstraints(2, y, GridBagConstraints.NORTHWEST));
			final JTextArea label2 = getLabelForOptionDescription(option);
			componentsPanel.add(label2, getGridBagConstraints(3, y, GridBagConstraints.NORTHWEST));
			y++;
		}
		// button to show the command
		final JButton buttonShowCommandLine = new JButton("Show command line with current parameters");
		buttonShowCommandLine.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				showCurrentCommandLine();

			}
		});
		buttonShowCommandLine.setToolTipText("Show command line with current parameters");
		final GridBagConstraints c2 = getGridBagConstraints(1, y, 1, 3, GridBagConstraints.WEST, 1.0);
		c2.fill = GridBagConstraints.NONE;
		componentsPanel.add(buttonShowCommandLine, c2);

		// now the run button
		final JButton button = new JButton("RUN");
		button.setToolTipText("Click to start program");
		final GridBagConstraints c = getGridBagConstraints(0, y, GridBagConstraints.CENTER);
		c.fill = GridBagConstraints.NONE;
		componentsPanel.add(button, c);

		// status at the south
		status = new JTextArea(10, 80) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 5204026945517631469L;

			@Override
			public void setText(String text) {
				super.setText(text);
				final String text2 = getText();
				if (text2 != null && !"".equals(text2)) {
					setCaretPosition(text2.length() - 1);
				}
			}

			@Override
			public void append(String text) {
				super.append(text);
				final String text2 = getText();
				if (text2 != null && !"".equals(text2)) {
					setCaretPosition(text2.length() - 1);
				}
			}
		};
		status.setEditable(false);
		status.setWrapStyleWord(true);
		status.setFont(new Font("Tahoma", Font.PLAIN, 11));
		status.setLineWrap(true);
		// set System.out to the textarea
		this.statusPrintStream = new MyPrintStream(new TextAreaOutputStream(status));

		//
		final JScrollPane scroll = new JScrollPane(status);
		scroll.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		final JPanel panelStatus = new JPanel(new BorderLayout());
		panelStatus.setBorder(BorderFactory.createTitledBorder("Status:"));
		panelStatus.add(scroll);
		panelStatus.setPreferredSize(new Dimension(400, 250));

		// split
		splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, scroll2, panelStatus);
		splitPane.setResizeWeight(0.5);
		splitPane.setOneTouchExpandable(true);
		getContentPane().add(splitPane, BorderLayout.CENTER);
		// what happens when pressing run
		button.addActionListener(getRunButtonAction(componentsByOption, program, status));

		// loadDefaults from defaults.properties
		loadDefaults(program.getCommandLineOptions());

		final Dimension preferredSize = new Dimension(SwingUtils.getFractionOfScreenWidthSize(0.8),
				Math.min(500, SwingUtils.getFractionOfScreenHeightSize(0.8)));
		setPreferredSize(preferredSize);
		SwingUtils.centerOnScreen(this);
		splitPane.setDividerLocation(0.8);
		//
		pack();
	}

	protected void showCurrentCommandLine() {
		String commandLineString;
		try {
			commandLineString = getCommandLineString(getCommandLineFromGui());
			showMessage("Command line from current options:");
			showMessage(commandLineString);
		} catch (final ParseException e) {
			e.printStackTrace();
			showError(e.getMessage());
		}

	}

	public MyPrintStream getStatusPrintStream() {
		return statusPrintStream;
	}

	protected String getCommandLineString(CommandLine commandLineFromGui) {
		final StringBuilder sb = new StringBuilder();
		final Option[] options = commandLineFromGui.getOptions();
		for (final Option option : options) {
			if (!"".equals(sb.toString())) {
				sb.append(" ");
			}
			sb.append("-" + option.getOpt());
			if (option.hasArg()) {
				sb.append(" ");
				final String value = option.getValue();
				if (value.contains(" ")) {
					sb.append("'");
				}
				sb.append(value);
				if (value.contains(" ")) {
					sb.append("'");
				}
			}
		}
		return sb.toString();
	}

	private void loadDefaults(Options options) {
		final AppDefaults defaults = new AppDefaults(componentsByOption.keySet());
		for (final String optionName : componentsByOption.keySet()) {
			final JComponent jComponent = componentsByOption.get(optionName);
			final String value = defaults.getPropertyValue(optionName);
			if (value != null) {
				if (jComponent instanceof JCheckBox) {
					((JCheckBox) jComponent).setSelected(Boolean.valueOf(value));
				} else if (jComponent instanceof JTextField) {
					((JTextField) jComponent).setText(value);
				}
			}
		}
	}

	private ActionListener getRunButtonAction(TMap<String, JComponent> componentsByOption,
			CommandLineProgramGuiEnclosable program, JTextArea status) {

		return new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {

				try {
					final long t1 = System.currentTimeMillis();
					final Options options = program.getCommandLineOptions();
					// log the options
					for (final String optionOpt : componentsByOption.keySet()) {
						final Option option = options.getOption(optionOpt);
						log.info(option.getOpt() + "=" + option.getValue());
					}
					// save defaults
					saveDefaults();
					clearStatus();
					final CommandLine commandLine = AutomaticGUICreator.this.getCommandLineFromGui();
					showMessage(getParametersString(commandLine));
					program.initTool(commandLine);

					showMessage("Parameters are correct. Starting program...");
					final ScheduledExecutorService service = new ScheduledThreadPoolExecutor(1);
					final Runnable command = new Runnable() {

						@Override
						public void run() {
							final ComponentEnableStateKeeper keeper = new ComponentEnableStateKeeper();
							keeper.addInvariableComponent(status);
							keeper.addInvariableComponent(status.getParent());
							keeper.addInvariableComponent(splitPane);
							try {
								keeper.keepEnableStates(AutomaticGUICreator.this);
								keeper.disable(AutomaticGUICreator.this);
								program.run();
								showMessage("Everything finished OK!!");
							} catch (final Exception e) {
								if (e.getCause() != null) {
									e.getCause().printStackTrace();
									showError(e.getCause().getMessage());
								} else {
									e.printStackTrace();
									showError(e.getMessage());
								}
							} finally {
								keeper.setToPreviousState(AutomaticGUICreator.this);
								final long t2 = System.currentTimeMillis();
								final long time = t2 - t1;
								showMessage("Process took " + DatesUtil.getDescriptiveTimeFromMillisecs(time));
							}
						}
					};
					service.schedule(command, 10, TimeUnit.MILLISECONDS);
				} catch (final SomeErrorInParametersOcurred e2) {
					showError(e2.getMessage());
				} catch (final Exception e1) {
					e1.printStackTrace();
					showError(e1.getMessage());
				}
			}

		};
	}

	public void showError(String message) {
		showMessage("Error: " + message);
		log.error(message);
	}

	public void showMessage(String message) {
		status.append(getFormattedTime() + ": " + message + "\n");
		log.info(message);
	}

	public static String getFormattedTime() {
		return format.format(new Date());

	}

	public static String getParametersString(CommandLine commandLine) {
		final StringBuilder sb = new StringBuilder("Program arguments:\n");
		final Option[] options = commandLine.getOptions();
		for (final Option option : options) {

			sb.append(option.getOpt());
			if (option.getValue() != null && !"".equals(option.getValue())) {
				sb.append(" = " + option.getValue());
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	protected void clearStatus() {
		status.setText("");
	}

	private void saveDefaults() {
		final AppDefaults appDefaults = new AppDefaults(componentsByOption.keySet());
		for (final String optionName : componentsByOption.keySet()) {
			final JComponent jComponent = componentsByOption.get(optionName);
			if (jComponent instanceof JCheckBox) {
				final Boolean checked = ((JCheckBox) jComponent).isSelected();
				appDefaults.setProperty(optionName, checked.toString());

			} else if (jComponent instanceof JTextField) {
				final String value = ((JTextField) jComponent).getText();
				appDefaults.setProperty(optionName, value);
			}
		}

	}

	public CommandLine getCommandLineFromGui() throws ParseException {
		final List<String> args = new ArrayList<String>();
		for (final String optionName : componentsByOption.keySet()) {

			final JComponent jComponent = componentsByOption.get(optionName);
			if (jComponent instanceof JCheckBox) {
				final boolean checked = ((JCheckBox) jComponent).isSelected();
				if (checked) {
					args.add("-" + optionName);
				}
			} else if (jComponent instanceof JTextField) {
				final String value = ((JTextField) jComponent).getText().trim();
				if (!"".equals(value)) {
					args.add("-" + optionName);
					args.add(value);
				}
			}
		}
		final CommandLineParser parser = new DefaultParser();
		final CommandLine cmd = parser.parse(program.getCommandLineOptions(), args.toArray(new String[0]));
		return cmd;
	}

	private static JComponent getComponentForOption(Option option, TMap<String, JComponent> componentsByOption) {
		JComponent ret = null;
		// if admits parameter use text box
		if (option.hasArg()) {
			final JTextField text = new JTextField();
			text.setColumns(30);
			ret = text;
		} else {
			// otherwise just a checkbox
			final JCheckBox check = new JCheckBox();
			ret = check;
		}
		componentsByOption.put(option.getOpt(), ret);
		ret.setToolTipText(option.getDescription());
		return ret;
	}

	private static String getNameFromOption(Option option) {
		String text = option.getOpt();
		if (option.hasLongOpt()) {
			text += " [" + option.getLongOpt() + "]";
		}
		return text;
	}

	private static JLabel getLabelForOptionName(Option option) {
		final String text = getNameFromOption(option);
		final JLabel label = new JLabel(text + ":");
		final Font font = label.getFont();
		label.setFont(new Font(font.getName(), Font.BOLD, font.getSize()));
		label.setToolTipText(option.getDescription());
		return label;
	}

	private static JTextArea getLabelForOptionDescription(Option option) {
		final String text = option.getDescription();
		final JTextArea label = new JTextArea(
//				"<html>" +
				text
//				+ "</html>"
		);
		label.setOpaque(false);
		label.setBorder(null);
		label.setEditable(false);
		label.setFocusable(false);
		label.setWrapStyleWord(true);
		label.setLineWrap(true);
		final Font font = new JLabel().getFont();
		label.setFont(new Font(font.getName(), Font.BOLD, font.getSize()));
		label.setToolTipText(option.getDescription());
		return label;
	}

	private static GridBagConstraints getGridBagConstraints(int x, int y, int anchor) {
		return getGridBagConstraints(x, y, 1, 1, anchor, 0.0);
	}

	private static GridBagConstraints getGridBagConstraints(int x, int y, int anchor, double weighty) {
		return getGridBagConstraints(x, y, 1, 1, anchor, weighty);
	}

	private static GridBagConstraints getGridBagConstraints(int x, int y, int gridheight, int gridwidth, int anchor,
			double weighty) {
		final GridBagConstraints gridBagConstraints = new GridBagConstraints();
		gridBagConstraints.gridx = x;
		gridBagConstraints.gridy = y;
		gridBagConstraints.insets = new Insets(5, 5, 5, 5);
		gridBagConstraints.gridheight = gridheight;
		gridBagConstraints.gridwidth = gridwidth;
		gridBagConstraints.weighty = weighty;
		gridBagConstraints.fill = GridBagConstraints.HORIZONTAL;
		gridBagConstraints.anchor = anchor;
		return gridBagConstraints;
	}

	public static AppVersion getVersion() {
		return getVersion(AppVersion.APP_PROPERTIES);
	}

	public static AppVersion getVersion(String propertiesFileName) {
		if (version == null) {
			try {
				final String tmp = PropertiesUtil
						.getProperties(new ClassPathResource(AppVersion.APP_PROPERTIES).getInputStream())
						.getProperty("assembly.dir");
				if (tmp.contains("v")) {
					version = new AppVersion(tmp.split("v")[1]);
				} else {
					version = new AppVersion(tmp);
				}
			} catch (final Exception e) {
//				e.printStackTrace();
			}
		}
		return version;

	}
}
