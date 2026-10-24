// Requests to the local grid-generator server.

// Decomposition searches and large exports can take a while, but a request
// that never answers must not leave the editor waiting forever.
const DEFAULT_TIMEOUT_MS = 120000;

export function getJson(url, options = {}) {
  return requestJson(url, { method: "GET" }, options);
}

export function postJson(url, body, options = {}) {
  return requestJson(
    url,
    {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        "X-SNaC-Grid-Request": "1",
      },
      body: JSON.stringify(body),
    },
    options
  );
}

// `signal` lets a caller cancel a request; a cancelled request rejects with the
// browser's AbortError so callers can recognize and ignore it.
async function requestJson(url, init, { allowInvalid = false, timeout = DEFAULT_TIMEOUT_MS, signal } = {}) {
  const controller = new AbortController();
  let timedOut = false;
  const timer = setTimeout(() => {
    timedOut = true;
    controller.abort();
  }, timeout);
  const cancel = () => controller.abort();
  signal?.addEventListener("abort", cancel, { once: true });
  let response;
  let text;
  try {
    response = await fetch(url, { ...init, signal: controller.signal });
    text = await response.text();
  } catch (error) {
    if (timedOut) throw new Error(`the server did not answer within ${Math.round(timeout / 1000)} s`);
    if (signal?.aborted) throw error;
    throw new Error("cannot reach the grid generator server; check that it is still running");
  } finally {
    clearTimeout(timer);
    signal?.removeEventListener("abort", cancel);
  }

  let payload = null;
  try {
    payload = JSON.parse(text);
  } catch (_error) {
    payload = null;
  }
  if (payload === null || typeof payload !== "object") {
    // Error pages from the server are HTML; show their message when present.
    const message = /<p>Message: ([^<]*?)\.?<\/p>/.exec(text ?? "")?.[1];
    if (!response.ok) throw new Error(`request failed (${response.status}${message ? `: ${message}` : ""})`);
    throw new Error("invalid server response");
  }
  if (!response.ok || (payload.ok === false && !allowInvalid)) {
    throw new Error(payload.error || payload.errors?.join("; ") || `request failed (${response.status})`);
  }
  return payload;
}
