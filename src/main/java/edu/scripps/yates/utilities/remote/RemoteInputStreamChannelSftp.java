package edu.scripps.yates.utilities.remote;

import java.io.IOException;
import java.io.InputStream;

import com.jcraft.jsch.ChannelSftp;
import com.jcraft.jsch.Session;
import com.jcraft.jsch.SftpException;

public class RemoteInputStreamChannelSftp extends InputStream {
	private final ChannelSftp sftp;
	private final String fileName;
	private Session session;
	private InputStream inputStream;

	public RemoteInputStreamChannelSftp(ChannelSftp sftp, Session session, String fileName) {
		this.sftp = sftp;
		this.fileName = fileName;
		this.session = session;
	}

	public InputStream getInputStream() throws SftpException {
		if (inputStream == null) {
			inputStream = sftp.get(fileName);
		}
		return inputStream;

	}

	@Override
	public void close() {
		sftp.exit();
		session.disconnect();
		try {
			getInputStream().close();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (SftpException e) {
			e.printStackTrace();
		}
	}

	@Override
	public int read() throws IOException {
		try {
			return getInputStream().read();
		} catch (SftpException e) {
			e.printStackTrace();
			throw new IOException(e);
		}
	}

}
