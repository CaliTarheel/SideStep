//! Minimal HTTP server for the globe viewer: `viewer <run_dir> [port]`.
//! Serves the embedded viewer page at `/`, the slice index at `/slices.json`, and run files under `/data/`.
use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::{Path, PathBuf};

const INDEX_HTML: &str = include_str!("../../viewer/index.html");

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let run = PathBuf::from(args.get(0).cloned().unwrap_or_else(|| "out/run".into()));
    let port: u16 = args.get(1).and_then(|p| p.parse().ok()).unwrap_or(8077);
    if !run.is_dir() { eprintln!("not a run directory: {}", run.display()); std::process::exit(1); }
    let listener = TcpListener::bind(("127.0.0.1", port)).expect("bind");
    eprintln!("viewer: http://127.0.0.1:{}/   (run: {})", port, run.display());
    for stream in listener.incoming() {
        if let Ok(stream) = stream {
            let run = run.clone();
            std::thread::spawn(move || handle(stream, &run));
        }
    }
}

fn handle(mut stream: TcpStream, run: &Path) {
    let mut buf = Vec::new();
    let mut tmp = [0u8; 4096];
    loop {
        match stream.read(&mut tmp) {
            Ok(0) => break,
            Ok(n) => { buf.extend_from_slice(&tmp[..n]); if buf.windows(4).any(|w| w == b"\r\n\r\n") || buf.len() > 65536 { break; } }
            Err(_) => return,
        }
    }
    let req = String::from_utf8_lossy(&buf);
    let path = req.lines().next().and_then(|l| l.split_whitespace().nth(1)).unwrap_or("/").to_string();
    let path = path.split('?').next().unwrap_or("/").to_string();

    if path == "/" || path == "/index.html" {
        respond(&mut stream, 200, "text/html; charset=utf-8", INDEX_HTML.as_bytes());
    } else if path == "/slices.json" {
        let mut times: Vec<i64> = std::fs::read_dir(run).map(|rd| rd.filter_map(|e| e.ok())
            .filter_map(|e| e.file_name().to_string_lossy().strip_prefix('t').and_then(|n| n.parse::<i64>().ok()))
            .collect()).unwrap_or_default();
        times.sort();
        let name = run.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
        let body = format!("{{\"run\":\"{}\",\"times\":[{}]}}", name, times.iter().map(|t| t.to_string()).collect::<Vec<_>>().join(","));
        respond(&mut stream, 200, "application/json", body.as_bytes());
    } else if let Some(rel) = path.strip_prefix("/data/") {
        if rel.contains("..") || rel.contains('\\') { respond(&mut stream, 403, "text/plain", b"forbidden"); return; }
        let file = run.join(rel);
        match std::fs::read(&file) {
            Ok(bytes) => {
                let ct = match file.extension().and_then(|e| e.to_str()).unwrap_or("") {
                    "png" => "image/png", "json" => "application/json", "csv" => "text/csv", _ => "application/octet-stream",
                };
                respond(&mut stream, 200, ct, &bytes);
            }
            Err(_) => respond(&mut stream, 404, "text/plain", b"not found"),
        }
    } else {
        respond(&mut stream, 404, "text/plain", b"not found");
    }
}

fn respond(stream: &mut TcpStream, code: u16, ctype: &str, body: &[u8]) {
    let status = match code { 200 => "OK", 403 => "Forbidden", 404 => "Not Found", _ => "Error" };
    let head = format!("HTTP/1.0 {} {}\r\nContent-Type: {}\r\nContent-Length: {}\r\nCache-Control: no-cache\r\nAccess-Control-Allow-Origin: *\r\nConnection: close\r\n\r\n", code, status, ctype, body.len());
    let _ = stream.write_all(head.as_bytes());
    let _ = stream.write_all(body);
    let _ = stream.flush();
}
