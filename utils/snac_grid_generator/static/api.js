export async function postJson(url, body, options = {}) {
  const response = await fetch(url, {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      "X-SNaC-Grid-Request": "1",
    },
    body: JSON.stringify(body),
  });
  const text = await response.text();
  let payload;
  try {
    payload = JSON.parse(text);
  } catch (_error) {
    throw new Error(response.ok ? "invalid server response" : `request failed (${response.status})`);
  }
  if (!response.ok || (payload.ok === false && !options.allowInvalid)) {
    throw new Error(payload.error || payload.errors?.join("; ") || "request failed");
  }
  return payload;
}
