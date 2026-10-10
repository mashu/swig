// Build the native aligner for this machine and stage it as native/swig, where
// cli-src/swig-cli-standalone.mjs embeds it (`npm run build:cli:standalone`).
// The release workflow stages each platform's binary from its own build instead.
import { execFileSync } from "node:child_process";
import { chmodSync, copyFileSync, mkdirSync } from "node:fs";
import { fileURLToPath } from "node:url";

const root = new URL("..", import.meta.url);
execFileSync("cargo", ["build", "--release", "--locked", "--manifest-path", fileURLToPath(new URL("engine/Cargo.toml", root)), "--bin", "swig"], { stdio: "inherit" });
const binary = new URL(`engine/target/release/${process.platform === "win32" ? "swig.exe" : "swig"}`, root);
const staged = new URL("native/swig", root);
mkdirSync(new URL("native/", root), { recursive: true });
copyFileSync(binary, staged);
chmodSync(staged, 0o755);
process.stdout.write(`Staged ${fileURLToPath(staged)}\n`);
