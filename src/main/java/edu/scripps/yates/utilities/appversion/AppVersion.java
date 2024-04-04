package edu.scripps.yates.utilities.appversion;

public class AppVersion implements Comparable<AppVersion> {
	private int mainVersion = -1;
	private int subVersion = -1;
	private int subsubVersion = -1;
	public static final String APP_PROPERTIES = "app.properties";

	public AppVersion(String versionString) {
		if (versionString == null || "".equals(versionString)) {
			throw new IllegalArgumentException("version string is null or empty");
		}
		final String version = getVersionString(versionString);
		if (version.contains(".")) {
			final String[] split = version.split("\\.");
			mainVersion = Integer.valueOf(split[0]);
			if (split.length > 1) {
				subVersion = Integer.valueOf(split[1]);
			}
			if (split.length > 2) {
				subsubVersion = Integer.valueOf(split[2]);
			}
		} else {
			mainVersion = Integer.valueOf(version);
		}
	}

	private String getVersionString(String versionString) {
		boolean isnumber = false;
		final StringBuilder ret = new StringBuilder();
		for (int i = 0; i < versionString.length(); i++) {
			final char c = versionString.charAt(i);
			if (Character.isDigit(c)) {
				isnumber = true;
				ret.append(c);
			} else if (c == '.' && isnumber) {
				ret.append(c);
			} else {
				isnumber = false;
			}
		}
		return ret.toString();
	}

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		if (mainVersion > -1) {
			sb.append(mainVersion);
		}
		if (subVersion > -1) {
			sb.append(".");
			sb.append(subVersion);
		}
		if (subsubVersion > -1) {
			sb.append(".");
			sb.append(subsubVersion);
		}
		return sb.toString();
	}

	@Override
	public int compareTo(AppVersion o) {
		if (o != null) {
			final int compare = Integer.compare(mainVersion, o.mainVersion);
			if (compare != 0) {
				return compare;
			} else {
				final int compare2 = Integer.compare(subVersion, o.subVersion);
				if (compare2 != 0) {
					return compare2;
				} else {
					return Integer.compare(subsubVersion, o.subsubVersion);

				}
			}
		}
		return 1;
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof AppVersion) {
			return ((AppVersion) o).compareTo(this) == 0;
		}
		return super.equals(o);
	}
}
