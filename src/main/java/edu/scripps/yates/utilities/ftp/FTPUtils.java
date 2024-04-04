package edu.scripps.yates.utilities.ftp;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.net.PrintCommandListener;
import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import org.apache.commons.net.ftp.FTPReply;
import org.apache.log4j.Logger;

import com.jcraft.jsch.Channel;
import com.jcraft.jsch.ChannelSftp;
import com.jcraft.jsch.ChannelSftp.LsEntry;
import com.jcraft.jsch.JSch;
import com.jcraft.jsch.JSchException;
import com.jcraft.jsch.Session;
import com.jcraft.jsch.SftpATTRS;
import com.jcraft.jsch.SftpException;
import com.jcraft.jsch.SftpProgressMonitor;

import edu.scripps.yates.utilities.files.FileUtils;

public class FTPUtils {
	private final static Logger log = Logger.getLogger(FTPUtils.class);
	private final static int timeoutInMillis = 10000;

	public static FTPClient loginFTPClient(String hostName, String userName, String password) throws IOException {
		return loginFTPClient(hostName, userName, password, -1);
	}

	public static FTPClient loginFTPClient(String hostName, String userName, String password, int port)
			throws IOException {

		final FTPClient ftpClient = new FTPClient();
		ftpClient.setDefaultTimeout(timeoutInMillis);
		// connect
		if (port > 0) {
			ftpClient.connect(hostName, port);
		} else {
			ftpClient.connect(hostName);
		}
		final int replyCode = ftpClient.getReplyCode();
		if (!FTPReply.isPositiveCompletion(replyCode)) {
			throw new IllegalArgumentException("Error while trying to connect to SFTP server at " + hostName
					+ ". Reply code: " + replyCode + " message:" + ftpClient.getReplyString());
		}
		showServerReply(ftpClient);
		// login
		final boolean success = ftpClient.login(userName, password);
		showServerReply(ftpClient);
		if (!success) {
			throw new IllegalArgumentException("Could not login to server " + hostName + " server");
		} else {
			log.debug("login succesful");
		}
		// enter in passive mode
		ftpClient.enterLocalPassiveMode();
		showServerReply(ftpClient);
		return ftpClient;

	}

	/**
	 * Download a file with regular FTP
	 * 
	 * @param host
	 * @param user
	 * @param pwd
	 * @param remoteFilePath
	 * @param outputStream
	 * @throws IOException
	 * @throws Exception
	 */
	public static void downloadFile(String host, String user, String pwd, String remoteFilePath,
			OutputStream outputStream) throws IOException {
		final FTPClient ftp = loginFTPClient(host, user, pwd);

		ftp.addProtocolCommandListener(new PrintCommandListener(new PrintWriter(System.out)));
		int reply;
		ftp.connect(host);
		reply = ftp.getReplyCode();
		if (!FTPReply.isPositiveCompletion(reply)) {
			ftp.disconnect();
			throw new IllegalArgumentException("Exception in connecting to FTP Server");
		}
		ftp.login(user, pwd);
		ftp.setFileType(FTP.BINARY_FILE_TYPE);
		ftp.enterLocalPassiveMode();
		ftp.retrieveFile(remoteFilePath, outputStream);
		if (ftp.isConnected()) {
			ftp.logout();
			ftp.disconnect();
		}
	}

	public static Session loginSSHClient(String hostName, String userName, String password, int port)
			throws JSchException {
		final JSch jsch = new JSch();
		final Session session = jsch.getSession(userName, hostName);
		session.setConfig("StrictHostKeyChecking", "no");
		// non-interactive version. Relies in host key being in known-hosts file
		session.setPassword(password);
		// connect
		session.connect();

		return session;

	}

	public static ChannelSftp openSFTPChannel(Session session) throws JSchException {
		final Channel channel = session.openChannel("sftp");
		channel.connect();
		final ChannelSftp sftpChannel = (ChannelSftp) channel;
		return sftpChannel;
	}

	public static void showServerReply(FTPClient ftpClient) {
		showServerReply(ftpClient, null);
	}

	public static void showServerReply(FTPClient ftpClient, PrintStream out) {
		final String[] replies = ftpClient.getReplyStrings();
		if (replies != null && replies.length > 0) {
			for (final String aReply : replies) {
				if (out != null) {
					out.println("SERVER: " + aReply);
				}
				log.debug("SERVER: " + aReply);
			}
		}
	}

	/**
	 * Creates a nested directory structure on a FTP server
	 * 
	 * @param ftpClient
	 *            an instance of org.apache.commons.net.ftp.FTPClient class.
	 * @param dirPath
	 *            Path of the directory, i.e /projects/java/ftp/demo
	 * @return true if the directory was created successfully, false otherwise
	 * @throws IOException
	 *             if any error occurred during client-server communication
	 */
	public static boolean makeDirectories(FTPClient ftpClient, String dirPath, PrintStream log) throws IOException {
		final String[] pathElements = dirPath.split("/");
		if (pathElements != null && pathElements.length > 0) {
			for (final String singleDir : pathElements) {
				if (!"".equals(singleDir)) {
					final boolean existed = ftpClient.changeWorkingDirectory(singleDir);
					if (!existed) {
						final boolean created = ftpClient.makeDirectory(singleDir);
						if (created) {
							log.println("CREATED directory: " + singleDir);
							ftpClient.changeWorkingDirectory(singleDir);
						} else {
							log.println("COULD NOT create directory: " + singleDir);
							showServerReply(ftpClient);
							throw new IllegalArgumentException("Error creating folder '" + singleDir + " ' in server");
						}
					}
				}
			}
		}
		return true;
	}

	public static long getSize(FTPClient ftpClient, String filePath) throws IOException {
		final FTPFile[] file = ftpClient.listFiles(filePath);
		if (file.length > 0) {
			return file[0].getSize();
		}
		return -1;
	}

	/**
	 * Determines whether a directory exists or not
	 * 
	 * @param dirPath
	 * @return true if exists, false otherwise
	 * @throws IOException
	 *             thrown if any I/O error occurred.
	 */
	public static boolean checkDirectoryExists(FTPClient ftpClient, String dirPath) throws IOException {
		ftpClient.changeWorkingDirectory(dirPath);
		final int returnCode = ftpClient.getReplyCode();
		if (returnCode == 550) {
			return false;
		}
		return true;
	}

	/**
	 * Determines whether a file exists or not
	 * 
	 * @param filePath
	 * @return true if exists, false otherwise
	 * @throws IOException
	 *             thrown if any I/O error occurred.
	 */
	public static boolean checkFileExists(FTPClient ftpClient, String filePath) throws IOException {
		final InputStream inputStream = ftpClient.retrieveFileStream(filePath);
		final int returnCode = ftpClient.getReplyCode();
		if (inputStream == null || returnCode == 550) {
			return false;
		}
		return true;
	}

	/**
	 * "Getting files recursively from remote host located under folder
	 * folderPath with a certain extension
	 * 
	 * @param extension
	 * @return
	 * @throws IOException
	 * @throws SftpException
	 */
	public static List<FTPFile> getFilesInFolderByExtension(FTPClient ftpClient, String folderPath, String extension)
			throws IOException {
		log.info("Getting files recursively from remote host located under folder: " + folderPath + " with extension "
				+ extension);
		final List<FTPFile> ret = new ArrayList<FTPFile>();
		final FTPFile[] listFiles = ftpClient.listFiles(folderPath);
		for (final FTPFile ftpFile : listFiles) {
			if (FilenameUtils.getExtension(ftpFile.getName()).equalsIgnoreCase(extension)) {
				ret.add(ftpFile);
			}
		}

		return ret;

	}

	public static List<String> getFileNamesInFolderByExtension(FTPClient ftpClient, String folderPath, String extension)
			throws IOException {
		final List<FTPFile> ftpFiles = getFilesInFolderByExtension(ftpClient, folderPath, extension);
		final List<String> ret = new ArrayList<String>();
		for (final FTPFile ftpFile : ftpFiles) {
			ret.add(ftpFile.getName());
		}
		return ret;
	}

	public static List<LsEntry> getFilesInFolderByExtension(Session sshClient, String folderPath, String extension)
			throws JSchException, SftpException {
		ChannelSftp sftpClient = null;
		try {
			sftpClient = openSFTPChannel(sshClient);
			final List<LsEntry> sftpFiles = sftpClient.ls(folderPath);
			final List<LsEntry> ret = new ArrayList<LsEntry>();
			for (final LsEntry sftpFile : sftpFiles) {
				if (FilenameUtils.getExtension(sftpFile.getFilename()).equalsIgnoreCase(extension)) {
					ret.add(sftpFile);
				}
			}
			return ret;
		} finally {
			if (sftpClient != null) {
				sftpClient.quit();
			}
		}
	}

	public static boolean exist(Session sshClient, String fileFullPath) throws JSchException, SftpException {
		ChannelSftp sftpClient = null;
		try {
			sftpClient = openSFTPChannel(sshClient);
			final List<LsEntry> sftpFiles = sftpClient.ls(fileFullPath);

			for (final LsEntry sftpFile : sftpFiles) {
				if (FilenameUtils.getName(sftpFile.getFilename())
						.equalsIgnoreCase(FilenameUtils.getName(fileFullPath))) {
					return true;
				}
			}
			return false;
		} finally {
			if (sftpClient != null) {
				sftpClient.quit();
			}
		}
	}

	public static long getSize(Session sshClient, String filePath) throws JSchException, SftpException {
		ChannelSftp sftpClient = null;
		try {
			sftpClient = openSFTPChannel(sshClient);
			final SftpATTRS stat = sftpClient.stat(filePath);
			if (stat != null) {
				return stat.getSize();
			}
		} finally {
			if (sftpClient != null) {
				sftpClient.quit();
			}
		}
		return -1;
	}

	public static long download(Session remoteServerSession, String fullPathToIP2, OutputStream outputStream,
			SftpProgressMonitor progressMonitor) throws JSchException, SftpException, IOException {
		long totalSize = 0;
		ChannelSftp sftpChannel = null;
		try {

			totalSize = FTPUtils.getSize(remoteServerSession, fullPathToIP2);

			if (totalSize > 0) {
				log.info("Downloading file '" + FilenameUtils.getName(fullPathToIP2) + " with size:"
						+ FileUtils.getDescriptiveSizeFromBytes(totalSize));
			}

			sftpChannel = FTPUtils.openSFTPChannel(remoteServerSession);

			sftpChannel.get(fullPathToIP2, outputStream, progressMonitor);

			outputStream.close();
			log.info("Transfer done.");
			return totalSize;

		} finally {
			if (sftpChannel != null) {
				sftpChannel.exit();
				sftpChannel.disconnect();
			} else {
				log.warn("Something happened!");
			}
		}

	}

	public static LsEntry getFileEntry(Session sshClient, String pathToFile) throws JSchException, SftpException {
		ChannelSftp sftpClient = null;
		try {
			sftpClient = openSFTPChannel(sshClient);
			final List<LsEntry> sftpFiles = sftpClient.ls(pathToFile);

			for (final LsEntry sftpFile : sftpFiles) {
				if (FilenameUtils.getName(sftpFile.getFilename()).equalsIgnoreCase(FilenameUtils.getName(pathToFile))) {
					return sftpFile;
				}
			}
		} finally {
			if (sftpClient != null) {
				sftpClient.quit();
			}
		}
		return null;
	}
}
