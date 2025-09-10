use std::env;
use std::process::Command;

fn main() {
    let version = env::var("CARGO_PKG_VERSION").unwrap();

    let git_output = Command::new("git")
        .args(["describe", "--tags", "--dirty", "--always"])
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| version.clone());

    println!("cargo:rustc-env=GIT_VERSION={}", git_output);
}
