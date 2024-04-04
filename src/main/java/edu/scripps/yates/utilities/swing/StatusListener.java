package edu.scripps.yates.utilities.swing;

public interface StatusListener<T> {

	public void onStatusUpdate(String statusMessage);

	public void onStatusUpdate(String statusMessage, T param);
}
