import http.server, os, sys, html
class H(http.server.SimpleHTTPRequestHandler):
    def list_directory(self, path):  # mimic nginx autoindex on ftp.ncbi.nlm.nih.gov
        entries = sorted(os.listdir(path)); rel = self.path
        out = [f"<html><head><title>Index of {rel}</title></head><body><h1>Index of {rel}</h1><hr><pre><a href=\"../\">../</a>"]
        for e in entries:
            full = os.path.join(path, e); isdir = os.path.isdir(full)
            name = e + ("/" if isdir else ""); size = "-" if isdir else str(os.path.getsize(full))
            out.append(f'<a href="{name}">{name}</a>{" " * max(1, 60-len(name))}25-Aug-2025 10:22 {size:>19}')
        out.append("</pre><hr></body></html>")
        body = "\n".join(out).encode(); self.send_response(200); self.send_header("Content-Type", "text/html"); self.send_header("Content-Length", str(len(body))); self.end_headers()
        import io; return io.BytesIO(body)
    def log_message(self, *a): pass
os.chdir(sys.argv[1]); http.server.ThreadingHTTPServer(("127.0.0.1", 8799), H).serve_forever()
