package edu.scripps.yates.utilities.jaxb.xpathquery;

public class XPathElement {
	public final String name;
	public final String filterValue;

	public XPathElement(String name, String filterValue) {
		super();
		this.name = name;
		if (!"".equals(filterValue)) {
			this.filterValue = filterValue;
		} else {
			this.filterValue = null;
		}
	}

	public String getName() {
		return name;
	}

	public String getFilterValue() {
		return filterValue;
	}

}
