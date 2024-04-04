package edu.scripps.yates.utilities.swing;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Rectangle;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JViewport;
import javax.swing.Scrollable;
import javax.swing.SwingConstants;

public class MyScrollablePanel extends JPanel implements Scrollable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -8751936348716442370L;
	private final JPanel content;

	public MyScrollablePanel(JPanel content) {
		super(new BorderLayout());
		this.content = content;
		add(content);
	}

	@Override
	public Dimension getPreferredScrollableViewportSize() {
		final Dimension preferredSize = content.getPreferredSize();
		if (getParent() instanceof JViewport) {
			preferredSize.width += ((JScrollPane) getParent().getParent()).getVerticalScrollBar()
					.getPreferredSize().width;
		}
		return preferredSize;
	}

	@Override
	public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) {
		return orientation == SwingConstants.HORIZONTAL ? Math.max(visibleRect.width / 10, 1)
				: Math.max(visibleRect.height / 10, 1);
	}

	@Override
	public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction) {
		return orientation == SwingConstants.HORIZONTAL ? Math.max(visibleRect.width * 9 / 10, 1)
				: Math.max(visibleRect.height * 9 / 10, 1);
	}

	@Override
	public boolean getScrollableTracksViewportWidth() {
		return true;
	}

	@Override
	public boolean getScrollableTracksViewportHeight() {
		if (getParent() instanceof JViewport) {
			final JViewport viewport = (JViewport) getParent();
			return getPreferredSize().height < viewport.getHeight();
		}
		return false;
	}
}
