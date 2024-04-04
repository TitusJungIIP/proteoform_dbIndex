package edu.scripps.yates.utilities.email;

import java.util.Properties;

import javax.mail.Message;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.AddressException;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;

import org.apache.log4j.Logger;

public class EmailSender {
	private static final Logger log = Logger.getLogger(EmailSender.class);

	public static String sendEmail(String subject, String body, String toemailAddress, String ccEmailAddress,
			String fromEmailAddress) {
		return sendEmail(subject, body, toemailAddress, ccEmailAddress, fromEmailAddress, false);
	}

	/**
	 * Send an email specifying the to and the from of the sender and receiver
	 * respectively, by using the localhost and mail.smtp.host as host
	 * 
	 * @param subject
	 * @param body
	 * @param toemailAddress
	 * @param fromEmailAddress
	 * @param formatAsHtml
	 * @return
	 */
	public static String sendEmail(String subject, String body, String toemailAddress, String ccEmailAddress,
			String fromEmailAddress, boolean formatAsHtml) {
		String errorMessage = null;
		// Recipient's email ID needs to be mentioned.
		final String to = toemailAddress;

		// Sender's email ID needs to be mentioned
		final String from = fromEmailAddress;

		// Assuming you are sending email from localhost
		final String host = "localhost";

		// Get system properties
		final Properties properties = System.getProperties();

		// Setup mail server
		properties.setProperty("mail.smtp.host", host);

		log.info("Trying to send an email from " + fromEmailAddress + " to " + toemailAddress + " with subject: "
				+ subject + " and body: " + body);
		// Get the default Session object.
		final Session session = Session.getDefaultInstance(properties);
		log.info("Email session created");
		try {
			// Create a default MimeMessage object.
			final Message message = new MimeMessage(session);

			// Set From: header field of the header.
			message.setFrom(new InternetAddress(from));

			// Set CC:
			if (ccEmailAddress != null) {
				message.setRecipient(Message.RecipientType.CC, new InternetAddress(ccEmailAddress));
			}

			// Set To: header field of the header.
			message.addRecipient(Message.RecipientType.TO, new InternetAddress(to));

			// Set Subject: header field
			message.setSubject(subject);

			// Now set the actual message
			if (formatAsHtml) {
				message.setContent(body, "text/html");
			} else {
				message.setText(body);
			}
			// Send message
			Transport.send(message);
			log.info("Email sent message successfully");
		} catch (final AddressException e) {
			errorMessage = "AddressException: " + e.getMessage();
		} catch (final javax.mail.MessagingException e) {
			errorMessage = "MessagingException: " + e.getMessage();
		}
		if (errorMessage != null)
			log.error("Email was not sent due to error: " + errorMessage);
		return errorMessage;
	}
}
