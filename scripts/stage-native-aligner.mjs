// Build the native aligner for this machine and stage it as native/swig, where
// cli-src/swig-cli-standalone.mjs embeds it (`npm run build:cli:standalone`).
// The release workflow stages each platform's binary from its own build instead.
import { execFileSync } from "node:child_process";
import { chmodSync, copyFileSync, mkdirSync } from "node:fs";
import { fileURLToPath } from "node:url";

const root = new URL("..", import.meta.url);
const manifest = fileURLToPath(new URL("engine/Cargo.toml", root));
// Cargo names the executable it built (CARGO_TARGET_DIR, build.target-dir and
// CARGO_BUILD_TARGET all move it), so stage that one and never a stale copy.
const messages = execFileSync(
  "cargo",
  ["build", "--release", "--locked", "--manifest-path", manifest, "--bin", "swig", "--message-format=json-render-diagnostics"],
  { encoding: "utf8", stdio: ["ignore", "pipe", "inherit"], maxBuffer: 256 * 1024 * 1024 },
);
const binary = messages
  .split("\n")
  .filter((line) => line.startsWith("{"))
  .map((line) => JSON.parse(line))
  .find((message) => message.reason === "compiler-artifact" && message.target?.name === "swig" && message.executable)?.executable;
if (!binary) throw new Error("cargo build did not report the swig executable");
const staged = new URL("native/swig", root);
mkdirSync(new URL("native/", root), { recursive: true });
copyFileSync(binary, staged);
chmodSync(staged, 0o755);
process.stdout.write(`Staged ${binary} as ${fileURLToPath(staged)}\n`);
