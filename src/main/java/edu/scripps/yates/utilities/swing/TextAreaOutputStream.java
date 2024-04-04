package edu.scripps.yates.utilities.swing;

import java.io.IOException;
import java.io.OutputStream;

import javax.swing.JTextArea;

/**
 * {@link JTextArea} that shows all the system out in the text.<br>
 * for this to work, you have to do something like:<br>
 * <code>JTextArea ta = new JTextArea(24, 80);<br>
                System.setOut(new PrintStream(new TextAreaOutputStream(ta)));
   </code>
 * 
 * @author salvador
 *
 */
public class TextAreaOutputStream extends OutputStream {

	private final JTextArea textArea;

	private final StringBuilder sb = new StringBuilder();

	public TextAreaOutputStream(final JTextArea textArea) {
		this.textArea = textArea;
	}

	@Override
	public void flush() {
	}

	@Override
	public void close() {
	}

	@Override
	public void write(int b) throws IOException {

		if (b == '\r') {
			return;
		}

		if (b == '\n') {
			final String text = sb.toString() + "\n";
			textArea.append(text);
			sb.setLength(0);
		} else {
			sb.append((char) b);
		}
	}

}
