use std::env;
use std::process::Command;

fn main() {
    let git_version = match env::var("SEQUINTOOLS_GIT_VERSION") {
        Ok(version) => version,
        Err(_) => {
            let version = env::var("CARGO_PKG_VERSION").unwrap();
            let git_output = Command::new("git")
                .args(["describe", "--tags", "--dirty", "--always"])
                .output()
                .ok()
                .and_then(|o| String::from_utf8(o.stdout).ok())
                .and_then(|s| {
                    let trimmed = s.trim();
                    if trimmed.is_empty() {
                        None
                    } else {
                        Some(trimmed.to_string())
                    }
                })
                .unwrap_or_else(|| version.clone());
            git_output
        }
    };
    println!("cargo:rustc-env=GIT_VERSION={}", git_version);
}
