package edu.scripps.yates.utilities.proteomicsmodel.factories;

import java.io.Serializable;
import java.util.Date;

import org.apache.commons.lang3.builder.HashCodeBuilder;

import edu.scripps.yates.utilities.proteomicsmodel.MSRun;
import edu.scripps.yates.utilities.proteomicsmodel.Project;

public class MSRunEx implements MSRun, Serializable {
	/**
	 *
	 */
	private static final long serialVersionUID = -247769229381526464L;
	private final String runID;
	private String path;
	private Date date;
	private int dbID;
	private Project project;
	private int hashCode = -1;

	public MSRunEx(String runID, String path) {
		if (runID == null)
			throw new IllegalArgumentException("Run ID is null");
		this.path = path;
		this.runID = runID;

	}

	protected String getRunID() {
		return runID;
	}

	/**
	 * @return the id of the run
	 */
	@Override
	public String getRunId() {
		return runID;
	}

	/**
	 * @return the name
	 */
	@Override
	public String getPath() {
		return path;
	}

	/**
	 * @return the date
	 */
	@Override
	public Date getDate() {
		return date;
	}

	/**
	 * @param date the date to set
	 */
	public void setDate(Date date) {
		this.date = date;
	}

	@Override
	public int getDBId() {
		return dbID;
	}

	public void setDBId(int id) {
		dbID = id;
	}

	@Override
	public Project getProject() {
		return project;
	}

	/**
	 * @param project the project to set
	 */
	public void setProject(Project project) {
		this.project = project;
	}

	@Override
	public void setPath(String runPath) {
		path = runPath;
	}

	@Override
	public String toString() {
		return runID;
	}

	@Override
	public int hashCode() {
		if (hashCode == -1) {
			hashCode = HashCodeBuilder.reflectionHashCode(runID);
		}
		return hashCode;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof MSRun) {
			final MSRun msRun = (MSRun) obj;
			return msRun.getRunId().equals(getRunId());
		}
		return super.equals(obj);
	}

}
