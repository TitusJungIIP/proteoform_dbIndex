package edu.scripps.yates.utilities.annotations.uniprot;

import org.apache.commons.lang.builder.HashCodeBuilder;

public class Gene {
	private final String name;
	private final String type;

	public Gene(String name, String type) {

		this.name = name;
		this.type = type;
	}

	public String getName() {
		return name;
	}

	public String getType() {
		return type;
	}

	@Override
	public String toString() {
		return "Gene [name=" + name + ", type=" + type + "]";
	}

	@Override
	public int hashCode() {
		return HashCodeBuilder.reflectionHashCode(toString());
	}

	@Override
	public boolean equals(Object obj) {
		return toString().equals(obj.toString());
	}

}
