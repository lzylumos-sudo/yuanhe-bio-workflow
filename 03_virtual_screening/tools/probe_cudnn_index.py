import ssl
import urllib.request

url = "https://developer.download.nvidia.com/compute/cudnn/redist/cudnn/linux-x86_64/"
request = urllib.request.Request(url, headers={"User-Agent": "zvss/0.9"})
with urllib.request.urlopen(request, timeout=60, context=ssl.create_default_context()) as response:
    body = response.read().decode("utf-8", "replace")
print("status: ok, bytes:", len(body))
print(body[:1200])
