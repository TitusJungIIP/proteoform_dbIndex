package edu.scripps.yates.utilities.swing;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.Window;

public class SwingUtils {
	public static int getFractionOfScreenHeightSize(double fraction) {
		final Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		final Double size = screenSize.getHeight() * fraction;
		return size.intValue();
	}

	public static int getFractionOfScreenWidthSize(double fraction) {
		final Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		final Double size = screenSize.getWidth() * fraction;
		return size.intValue();
	}

	public static Dimension getScreenDimension() {
		return Toolkit.getDefaultToolkit().getScreenSize();
	}

	public static void setComponentPreferredSizeRelativeToScreen(Component component, double widthFraction,
			double heightFraction) {
		final Dimension preferredSize = new Dimension(SwingUtils.getFractionOfScreenWidthSize(widthFraction),
				SwingUtils.getFractionOfScreenHeightSize(heightFraction));
		component.setPreferredSize(preferredSize);
	}

	public static void centerOnScreen(Window window) {
		window.pack();
		final java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
		final java.awt.Dimension dialogSize = window.getSize();
		window.setLocation((screenSize.width - dialogSize.width) / 2, (screenSize.height - dialogSize.height) / 2);
	}
}
