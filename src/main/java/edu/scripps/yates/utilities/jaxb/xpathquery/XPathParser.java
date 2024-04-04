package edu.scripps.yates.utilities.jaxb.xpathquery;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class XPathParser implements Iterator<XPathElement>, Cloneable {
	private final List<XPathElement> elements = new ArrayList<XPathElement>();

	private XPathElement currentElement;
	private final String xpath;

	private Iterator<XPathElement> iterator;

	public XPathParser(String xpath) {
		this.xpath = xpath;
		parse(xpath);
	}

	private void parse(String xpath) {

		final List<String> elements2 = new ArrayList<String>();
		if (xpath.contains("/")) {
			final String[] split = xpath.split("/");
			for (final String element : split) {
				elements2.add(element.trim());
			}
		} else {
			elements2.add(xpath.trim());
		}

		for (final String element : elements2) {
			if (element.contains("$")) {
				final String[] split = element.split("\\$");
				elements.add(new XPathElement(split[0].trim(), null));

				final String attribute = split[1].trim();
				if (attribute != null) {
					if (attribute.contains("=")) {
						final String[] split2 = attribute.split("=");
						elements.add(new XPathElement(split2[0].trim(), split2[1].trim()));

					}
				}

			} else {
				elements.add(new XPathElement(element.trim(), null));

			}
		}
		iterator = elements.iterator();
	}

	public List<XPathElement> getElements() {
		return elements;
	}

	@Override
	public boolean hasNext() {
		return iterator.hasNext();
	}

	@Override
	public XPathElement next() {
		final XPathElement next = iterator.next();
		currentElement = next;
		return next;
	}

	public XPathElement getCurrentElement() {
		return currentElement;
	}

	public String getFilterValue() {
		return currentElement.getFilterValue();
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		final XPathParser ret = new XPathParser(xpath);
		return ret;
	}
}
