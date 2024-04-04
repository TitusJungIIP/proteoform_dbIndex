package edu.scripps.yates.utilities.http;

import java.io.IOException;

import org.apache.commons.io.IOUtils;
import org.apache.http.HttpEntity;
import org.apache.http.HttpHost;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.client.methods.HttpRequestBase;
import org.apache.http.client.protocol.HttpClientContext;
import org.apache.http.conn.routing.HttpRoute;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.impl.conn.PoolingHttpClientConnectionManager;
import org.apache.log4j.Logger;

public class ThreadSafeHttpClient {
	private final static Logger log = Logger.getLogger(ThreadSafeHttpClient.class);
	private static PoolingHttpClientConnectionManager cm;

	private static PoolingHttpClientConnectionManager getConnectionManager() {
		if (cm == null) {
			cm = new PoolingHttpClientConnectionManager();
			// Increase max total connection to 200
			cm.setMaxTotal(200);
			// Increase default max connection per route to 20
			cm.setDefaultMaxPerRoute(20);
			// Increase max connections for localhost:80 to 50
			final HttpHost localhost = new HttpHost("locahost", 80);
			cm.setMaxPerRoute(new HttpRoute(localhost), 50);
		}
		return cm;
	}

	public static CloseableHttpClient createNewHttpClient() {
		return HttpClients.custom().setConnectionManager(getConnectionManager()).build();

	}

	/**
	 * 
	 * @param httpClient     can be created with <code>createNewHttpClient()</code>
	 * @param httpUriRequest it can be a {@link HttpGet} or {@link HttpPost}
	 * @return
	 * @throws ClientProtocolException
	 * @throws IOException
	 */
	public static String getStringResponse(CloseableHttpClient httpClient, HttpRequestBase httpUriRequest)
			throws ClientProtocolException, IOException {
		try {
			final CloseableHttpResponse response = httpClient.execute(httpUriRequest, HttpClientContext.create());
			if (response.getStatusLine().getStatusCode() == 200) {
				final HttpEntity entity = response.getEntity();
				final String ret = IOUtils.toString(entity.getContent());
				return ret;
			}
			log.error("Error from http client. Status code:" + response.getStatusLine().getStatusCode() + ", Reason:"
					+ response.getStatusLine().getReasonPhrase());
			return null;
		} finally {
			httpUriRequest.releaseConnection();
		}

	}

	/**
	 * shutdown the singleton connection manager and set it to null so that next
	 * time a new HTTPClient is created, a new connection manager will be
	 * created.<br>
	 * call this method after closing your client.
	 */
	public static void closeConnectionManager() {
		getConnectionManager().shutdown();
		cm = null;
	}
}
