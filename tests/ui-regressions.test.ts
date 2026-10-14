import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

import { logoGlyphRun } from "../src/logo-glyphs.ts";
import { AMINO_ACID_LOGO_SYMBOLS, NUCLEOTIDE_LOGO_SYMBOLS } from "../src/probability-logo.ts";

test("clickable individual-sequence views retain CDR3 context", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const post = fs.readFileSync(new URL("../src/post-analysis.tsx", import.meta.url), "utf8");

  // Paged AIRR records show AA with nucleotide fallback and retain nucleotide
  // context when both forms are present.
  assert.match(app, /row\.cdr3Aa \|\| row\.cdr3 \|\| "—"/);
  assert.match(app, /row\.cdr3Aa && row\.cdr3/);

  // The sparse VDDJ explorer hydrates CDR3 from the main index rather than
  // expecting it to be duplicated in the double-D sidecar.
  assert.match(app, /cdr3:record\.values\.cdr3\|\|index\?\.cdr3/);
  assert.match(app, /<th>CDR3<\/th>/);

  // Chimera candidates and lineage-member click targets also expose CDR3.
  assert.match(post, /indexed\?\.cdr3Aa\|\|indexed\?\.cdr3\|\|"CDR3 —"/);
  assert.match(post, /row\.values\.cdr3_aa\|\|row\.values\.cdr3\|\|"CDR3 —"/);
});

test("FASTQ quality control is an initially collapsed pre-assignment step", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");

  // No `open` property means the native details box is initially collapsed.
  assert.match(app, /<details className=\{`fastq-quality-control/);
  assert.doesNotMatch(app, /<details[^>]*fastq-quality-control[^>]*\sopen(?:=|\s|>)/);
  assert.match(app, /3′ trim → minimum retained length → expected-error filter → optional random subsample → V\(D\)J assignment/);
  assert.match(app, /FASTA and AIRR records will pass through this step unchanged/);
});

test("the web entry point exposes three focused workflows and keeps UCA trees observed-only", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const post = fs.readFileSync(new URL("../src/post-analysis.tsx", import.meta.url), "utf8");
  const uca = fs.readFileSync(new URL("../src/phylo-uca/panel.tsx", import.meta.url), "utf8");
  const styles = fs.readFileSync(new URL("../src/globals.css", import.meta.url), "utf8");

  assert.match(app, /VDJ assignment, annotation, and visualization/);
  assert.match(app, /End-to-end repertoire analysis/);
  assert.match(app, /Single lineage analysis/);
  assert.match(app, /Swig runs on your machine\. Your sequence data are not uploaded/);
  assert.doesNotMatch(app, /Local annotation of/);
  assert.match(app, /useState<CallingProfile>\("r_optimized"\)/);
  assert.doesNotMatch(app, /AssignerStrategy|assignerStrategy/);
  assert.match(app, /!focusedWebMode && <section className=\{`analysis-card pipeline-card/);
  assert.match(app, /directLineage=\{session\.webMode === "lineage"\}/);
  assert.match(styles, /\.workflow-choice-grid[\s\S]*grid-template-columns: repeat\(3/);

  assert.match(post, /Upload Newick/);
  assert.match(post, /allowGermlineGuide|prepareUploadedLineageTree\(text, observedNames, true\)/);
  assert.match(post, /treeRun\?\.origin === "uploaded" && treeRun\.observedOnly/);
  assert.match(post, /Contains germline guide · display only; UCA will infer a fresh observed-only tree/);
  assert.match(uca, /Use uploaded observed-only tree/);
  assert.match(uca, /Infer a fresh tree with FastTree/);
  assert.match(uca, /prepareObservedOnlyAlignment\(alignment, GERMLINE_OUTGROUP\)/);
  assert.match(uca, /usingUploadedTree \? suppliedObservedTree : inferredObservedTree!\.newick/);
});

test("concatenated gzip uploads require an explicit merged-versus-separate sample choice", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const gzip = fs.readFileSync(new URL("../src/gzip-members.ts", import.meta.url), "utf8");

  assert.match(app, /CONCATENATED GZIP · SAMPLE STRUCTURE/);
  assert.match(app, /Import as \{pendingGzipImport\.members\.length\} separate samples/);
  assert.match(app, /Merge into one sample/);
  assert.match(app, /Independent sample, donor, cohort, timepoint, and tissue fields/);
  assert.match(app, /gzipMemberSource\(candidate\.file,members\)/);
  assert.match(gzip, /startsAtRecordBoundary/);
  assert.match(gzip, /multi-member candidates are each decompressed once/);
});

test("session restoration applies saved metadata during the first index pass", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const restore = app.slice(app.indexOf("async function restoreSavedSession"), app.indexOf("async function restoreLinkedAirr"));
  const configure = restore.indexOf("configureStudyMetadataForImport(restoredDatasets)");
  const append = restore.indexOf("store.appendBatch(header,body)");
  assert.ok(configure >= 0 && configure < append, "saved metadata must be installed before AIRR indexing starts");
  assert.doesNotMatch(restore, /updateStudyMetadata/);
  assert.doesNotMatch(restore, /Applying saved study metadata to local indexes/);
});

test("analysis-scale running states expose cancellation and every interface family has a method document", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const post = fs.readFileSync(new URL("../src/post-analysis.tsx", import.meta.url), "utf8");
  const uca = fs.readFileSync(new URL("../src/phylo-uca/panel.tsx", import.meta.url), "utf8");
  assert.match(app, /Cancel analysis/);
  assert.match(app, /Cancel reference preparation/);
  assert.match(app, /Cancel inspection/);
  assert.match(app, /Cancel re-indexing/);
  assert.match(app, /Cancel results export/);
  assert.match(app, /Cancel session save/);
  assert.match(app, /Cancel search/);
  assert.match(app, /cancelSessionRestore/);
  assert.match(post, /runtime\.beginTransaction\(\)/);
  assert.match(post, /runtime\.cancelAndRestore/);
  assert.match(post, /Restoring…/);
  assert.match(uca, /abortRef\.current\?\.abort\(\)/);

  const index = fs.readFileSync(new URL("../public/METHODS_INDEX.md", import.meta.url), "utf8");
  for (const document of [
    "01_INPUT_QC_AND_STUDY.md", "02_REFERENCE_PREPARATION.md", "03_VDJ_ASSIGNMENT.md",
    "04_DOUBLE_D_SCREEN.md", "05_STORAGE_DASHBOARD_SELECTION.md", "06_COLLAPSE_AND_DENOISING.md",
    "07_CHIMERA_INFERENCE.md", "08_LINEAGE_ASSIGNMENT.md", "09_SHM_AND_GERMLINE_DIAGNOSTICS.md",
    "10_RETRIEVAL_AND_NEIGHBOURS.md", "11_ALIGNMENT_AND_PHYLOGENY.md",
    "12_EXECUTION_PERSISTENCE_CANCELLATION.md",
  ]) {
    assert.match(index, new RegExp(document.replaceAll(".", "\\.")));
    assert.equal(fs.existsSync(new URL(`../public/methods/${document}`, import.meta.url)), true, `missing ${document}`);
  }
  assert.match(index, /REPERTOIRE_ALLELE_REFINEMENT\.md/);
  assert.match(index, /PHYLO_UCA_INFERENCE\.md/);
});

test("personalized germline inference is isolated beside missing-V diagnostics and selects lowest current V-SHM",()=>{
  const post=fs.readFileSync(new URL("../src/post-analysis.tsx",import.meta.url),"utf8");
  const engine=fs.readFileSync(new URL("../src/personalized-germline.ts",import.meta.url),"utf8");
  assert.match(post,/Infer expressed V allele set/);
  assert.match(post,/Lowest current V-SHM/);
  assert.match(post,/no circularity correction/);
  assert.match(engine,/observation\.shmRate < existing\.shmRate/);
  assert.doesNotMatch(engine,/random|softmax|gumbel/i);
});

test("assignment progress uses stable geometry and live scheduler telemetry", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const scheduler = fs.readFileSync(new URL("../src/aligner-worker.ts", import.meta.url), "utf8");
  const runtime = fs.readFileSync(new URL("../src/aligner-runtime.ts", import.meta.url), "utf8");
  const styles = fs.readFileSync(new URL("../src/globals.css", import.meta.url), "utf8");

  assert.match(app, /Live assignment telemetry/);
  assert.match(app, /Queries \/ s/);
  assert.match(app, /Active workers/);
  assert.match(app, /Draining final batches/);
  assert.match(app, /onTelemetry: setAssignmentTelemetry/);
  assert.match(runtime, /onTelemetry\?: \(telemetry: AssignmentTelemetry\)/);
  assert.match(scheduler, /activeInitializers -= 1/);
  assert.match(scheduler, /postInitializationTelemetry\(request\.id, workerCount, activeInitializers\)/);
  assert.match(scheduler, /activeWorkers: slots\.reduce/);
  assert.match(scheduler, /acknowledged \+= result\.count/);
  assert.match(scheduler, /ASSIGNMENT_TELEMETRY_INTERVAL_MS/);
  assert.match(styles, /\.progress-orbit > span b[\s\S]*width: 3ch/);
  assert.match(styles, /\.progress-copy h2[\s\S]*height: 3\.45em/);
  assert.match(styles, /\.progress-copy h2[\s\S]*font-variant-numeric: tabular-nums/);
  assert.match(styles, /\.assignment-telemetry[\s\S]*grid-template-columns: repeat\(3/);
  assert.match(styles, /\.assignment-telemetry > div[\s\S]*grid-template-rows: 14px 34px 24px/);
  assert.match(styles, /\.assignment-telemetry dd[\s\S]*font-variant-numeric: tabular-nums/);
});

test("assignment is one progressive action page while independent result tools remain contextual", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const repertoire = fs.readFileSync(new URL("../src/repertoire-charts.tsx", import.meta.url), "utf8");
  const post = fs.readFileSync(new URL("../src/post-analysis.tsx", import.meta.url), "utf8");
  const css = fs.readFileSync(new URL("../src/globals.css", import.meta.url), "utf8");
  const cli = fs.readFileSync(new URL("../cli-src/swig-cli.mjs", import.meta.url), "utf8");
  const cliWorker = fs.readFileSync(new URL("../cli-src/swig-worker.mjs", import.meta.url), "utf8");
  const computeWorker = fs.readFileSync(new URL("../src/aligner-compute-worker.ts", import.meta.url), "utf8");
  const aligner = fs.readFileSync(new URL("../src/aligner-engine.ts", import.meta.url), "utf8");

  assert.match(app, /analysis-layout single-action-layout/);
  assert.doesNotMatch(app, /AnalysisWorkspace|analysisWorkspace|Analysis setup sections/);
  assert.doesNotMatch(app, /settings-strip/);
  assert.match(app, /<b>Advanced options<\/b>/);
  assert.match(app, /<option value="r_optimized">R-optimized · default<\/option>/);
  assert.match(app, /<option value="sensitive_d">Sensitive-D<\/option>/);
  assert.match(app, /one placement DAG/);
  assert.match(app, /V \+2\/−3\/−9\/−1, D \+2\/−3\/−13\/−1 with a 4-nt exact-run floor/);
  assert.match(cli, /vdjTuning\(options\)/);
  // Browser and command-line workers share one wrapper of the Rust module.
  assert.match(cliWorker, /import \{ AlignerEngine \} from "\.\.\/src\/aligner-engine\.ts"/);
  assert.match(computeWorker, /AlignerEngine\.create\(/);
  assert.match(aligner, /\["profile", settings\.callingProfile\]/);
  assert.match(aligner, /WebAssembly\.instantiate\(module, \{\}\)/);
  assert.match(app, /Customize individual loci, V\/D\/J\/C sources, or allele inclusion/);
  assert.match(app, /aria-label="Sequence result panels"/);
  assert.match(app, /setSequenceWorkspace\("detail"\)/);
  assert.match(app, /Records \+ filters/);
  assert.match(app, /sequence-records-stack/);
  assert.doesNotMatch(app, /useState<"filters"\|"records"/);
  assert.match(app, /context-rail-secondary/);
  assert.match(app, /<span>Double-D explorer<small>/);
  assert.doesNotMatch(app, /results-tab-double-d|results-panel-double-d/);
  assert.match(app, /results-hero" hidden aria-hidden="true"/);
  assert.match(app, /className="results-rail-tools"/);
  assert.match(app, /sidebarTools=\{resultsSidebarTools\(\)\}/);
  assert.match(app, /overview=\{repertoireOverview\}/);
  assert.match(app, /className="results-tool-drawer"/);
  assert.match(app, /<strong>Study design<\/strong>/);
  assert.match(app, /<strong>Sample colors<\/strong>/);
  assert.match(app, /<b>01<\/b><span>Repertoire/);
  assert.match(app, /FacetPicker label=\{`\$\{segment\} gene or allele`\}/);
  assert.match(app, /Include multi-call assignments containing this/);
  assert.match(app, /addFieldHelp\(root\)/);
  assert.match(app, /<select value=\{species\?\.name\?\?""\}/);
  assert.doesNotMatch(app, /speciesDraft/);
  assert.match(app, /setPendingDatabaseId\(nextId\)/);
  assert.match(app, /aria-busy=\{busy\} value=\{value\} onChange=/);
  assert.match(app, /Covered C-prefix identity gate/);
  assert.match(app, /A failed C gate clears only C fields; V\/D\/J calls and boundaries are retained/);
  assert.match(app, /<option value="gzip">Gzip \(\.gz\)<\/option>/);

  assert.match(repertoire, /aria-label="Repertoire panels"/);
  assert.match(repertoire, /panel==="usage"/);
  assert.doesNotMatch(repertoire, /panel==="controls"|setPanel\("controls"\)/);
  assert.match(repertoire, /repertoire-global-settings/);
  assert.match(repertoire, /\{sidebarTools\}/);
  assert.match(repertoire, /\{overview\}/);

  assert.match(post, /type PostWorkspaceId = "overview" \| PostModuleId/);
  assert.match(post, /aria-label="Post-analysis sections"/);
  assert.match(post, /new Set<PostModuleId>\(activeWorkspace === "overview" \? \[\] : \[activeWorkspace\]\)/);
  assert.match(post, /Advanced \{collapseMode === "fad"/);
  assert.match(post, /Advanced call matching and performance settings/);
  assert.match(post, /Advanced missing-allele evidence thresholds/);
  assert.match(post, /Study metadata filters/);
  assert.match(post, /Additional D, constant-gene, isotype, and rare-event filters/);
  assert.match(post, /FacetPicker label="Sample"/);
  assert.match(post, /FacetPicker label="V gene or allele"/);
  assert.match(post, /vCallIncludeAmbiguous/);
  assert.match(post, /<span>Allele pooling<small>/);
  assert.match(post, /Resolve ambiguous germline calls by pooling repertoire evidence/);
  assert.match(post, /<AlleleRefinementPanel references=\{references\} options=/);
  assert.match(post, /bounded native compression stream/);
  assert.match(post, /\{sidebarTools\}/);
  assert.match(fs.readFileSync(new URL("../src/allele-refinement/panel.tsx", import.meta.url), "utf8"), /Advanced evidence-kernel and model settings/);
  assert.match(css, /\.post-context-main > \.post-module\.is-collapsed \{ display: none !important; \}/);
  assert.match(css, /\.pipeline-stage-grid article:not\(\.enabled\) \.pipeline-fields/);
  assert.match(css, /\.post-module > label\.constant-collapse-policy/);
  assert.match(css, /\.results-table col\.column-cdr3/);
  assert.match(css, /\.sequence-context-main \.detail-shell[\s\S]*overflow-x: auto/);
  assert.match(css, /results are an application workspace, not a report masthead/);
  assert.match(css, /\.results-application-page \.results-view-tabs[\s\S]*background: var\(--ink\)/);
  assert.match(css, /\.post-analysis-heading \{ display: none !important; \}/);
});

test("allele pooling exposes parameter-responsive reference and hard assignment diagnostics", () => {
  const panel = fs.readFileSync(new URL("../src/allele-refinement/panel.tsx", import.meta.url), "utf8");
  const views = fs.readFileSync(new URL("../src/allele-refinement/diagnostic-views.tsx", import.meta.url), "utf8");
  const diagnostics = fs.readFileSync(new URL("../src/allele-refinement/diagnostics.ts", import.meta.url), "utf8");
  const post = fs.readFileSync(new URL("../src/post-analysis.tsx", import.meta.url), "utf8");
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  assert.match(panel, /ReferenceKernelInspector references=\{references\} options=\{options\}/);
  assert.match(panel, /AlleleAssignmentShiftChart results=\{segmentResults\}/);
  assert.match(panel, /Best posterior if confidence passes/);
  assert.match(panel, /Best posterior for every modeled record/);
  assert.match(panel, /Fast hurdle active set/);
  assert.match(panel, /Continuous Dirichlet mixture/);
  assert.match(panel, /Inclusion posterior threshold/);
  assert.doesNotMatch(panel, /allele-model-table/);
  assert.doesNotMatch(panel, /model\.alleles\.slice\(0,8\)/);
  assert.match(views, /label="Reference allele"/);
  assert.match(views, />Differences<\/button>/);
  assert.match(views, /row\.primary \? <small>bar omitted<\/small>/);
  assert.match(views, /Best-match allele counts before and after reassignment/);
  assert.match(views, /Vanished alleles only/);
  assert.match(views, /local_best_count/);
  assert.match(views, /SVG ↓/);
  assert.match(views, /Download surviving allele reference/);
  assert.match(views, /minimumReferenceReads/);
  assert.match(diagnostics, /unstripped\.some\(\(row\) => row\.sequence\[column\] !== "-"\)/);
  assert.match(diagnostics, /export function hardAssignmentShiftData/);
  assert.match(diagnostics, /export function survivingAlleleReference/);
  assert.match(diagnostics, /vanishes: beforeCount > 0 && afterCount === 0/);
  assert.match(panel, /Cross-donor override/);
  assert.match(panel, /Evidence never crosses participants/);
  assert.match(post, /alleleRuntime\.run\(store, references, alleleOptions, null,/);
  const automaticAlleles = post.indexOf("if (autoPipeline.alleleRefinement.enabled)");
  const automaticCollapse = post.indexOf("if (autoPipeline.collapse.enabled)");
  assert.ok(automaticAlleles >= 0 && automaticAlleles < automaticCollapse, "automatic pooling must run before collapse");
  const manualCollapse = post.slice(post.indexOf("async function runDedup()"), post.indexOf("async function applyDedupFilter()"));
  assert.doesNotMatch(manualCollapse, /setRepertoireCallOverrides\(null\)/);
  assert.match(app, /1 · Repertoire allele pooling/);
  assert.match(app, /Uses policy-selected V\/D\/J calls/);
});

test("post-analysis skipping, richer lineage rows, CLI export, and lazy lineage studies remain exposed",()=>{
  const post=fs.readFileSync(new URL("../src/post-analysis.tsx",import.meta.url),"utf8");
  const postWorker=fs.readFileSync(new URL("../src/post-analysis-worker.ts",import.meta.url),"utf8");
  const partitionWorker=fs.readFileSync(new URL("../src/post-analysis-partition-worker.ts",import.meta.url),"utf8");
  const app=fs.readFileSync(new URL("../src/swig-app.tsx",import.meta.url),"utf8");
  const lineageStudy=fs.readFileSync(new URL("../src/lineage-study-page.tsx",import.meta.url),"utf8");
  const standalone=fs.readFileSync(new URL("../cli-src/swig-cli-standalone.mjs",import.meta.url),"utf8");
  const release=fs.readFileSync(new URL("../.github/workflows/release-cli.yml",import.meta.url),"utf8");
  assert.match(post,/autoPipeline\?\.selection\.enabled\?\[\]:\["selection"\]/);
  assert.match(post,/Skip step/);
  assert.match(post,/Include step/);
  assert.match(post,/Contains reads from sample/);
  assert.match(post,/Other samples may also occur in the lineage/);
  assert.match(post,/sample\.sampleId\}: \{sample\.abundance\.toLocaleString\(\)\} read/);
  assert.match(post,/Mean SHM/);
  assert.match(post,/SHM upper q95/);
  assert.match(post,/Representative CDR3/);
  assert.match(post,/Open selected together/);
  assert.match(post,/Open together/);
  assert.match(post,/Indel-aware edit distance/);
  assert.match(post,/Require compatible V or J/);
  assert.match(post,/Post-lineage analysis/);
  assert.match(postWorker,/new Worker\(new URL\("\.\/post-analysis-partition-worker\.ts"/);
  assert.match(postWorker,/records\.length >= 10_000/);
  assert.match(postWorker,/>= 500 \? request\.workers : 1/);
  assert.match(partitionWorker,/runDenoisePartitionJob/);
  assert.match(partitionWorker,/runExactDedupJob/);
  assert.match(app,/Export CLI config/);
  assert.match(app,/Run this setup with swig-cli/);
  assert.match(app,/exportConfiguredCliConfig/);
  assert.match(app,/Cancel config export/);
  assert.match(app,/Cancel CLI config export/);
  assert.match(app,/Load lineage study/);
  assert.match(lineageStudy,/readLineageAirrSlice/);
  assert.match(lineageStudy,/Promise\.all\(selected\.map/);
  assert.match(lineageStudy,/Open together/);
  assert.match(lineageStudy,/Representative CDR3/);
  assert.match(lineageStudy,/Samples and read counts/);
  assert.match(standalone,/with \{ type: "file" \}/);
  assert.match(standalone,/runCli\(\{wasmPath,referencePackPath,nativeAlignerPath\}\)/);
  assert.match(standalone,/import nativeAlignerPath from "\.\.\/native\/swig" with \{ type: "file" \}/);
  assert.match(release,/needs: \[wasm, aligner\]/);
  assert.match(release,/aligner: swig-aligner-windows-x64\.exe/);
  // Bun 1.4.2 closes node:readline under a paused consumer ("readline was closed"); stay on 1.3.14.
  assert.match(release,/bun-version: 1\.3\.14/);
  assert.match(release,/Build release WebAssembly core/);
  assert.match(release,/release-swig-wasm/);
  assert.match(release,/build-swig-wasm\.sh/);
  assert.match(release,/wasm32-unknown-unknown/);
  for(const target of ["bun-linux-x64-baseline","bun-linux-arm64","bun-darwin-x64","bun-darwin-arm64","bun-windows-x64-baseline"]){assert.match(release,new RegExp(target));}
  assert.match(release,/Smoke-test the actual executable/);
  // Linux aligners link glibc statically (swig-cli runs from glibc 2.17); Windows runs its executable.
  assert.match(release,/target-feature=\+crt-static/);
  assert.match(release,/needs: \[build, aligner, windows-smoke\]/);
  const cargoConfig=fs.readFileSync(new URL("../.cargo/config.toml",import.meta.url),"utf8");
  assert.match(cargoConfig,/\[target\.x86_64-pc-windows-msvc\]\nrustflags = \["-C", "target-feature=\+crt-static"\]/);
});

test("lineage alignment offers direct Alivibe-compatible MSA and retains manual round trip", () => {
  const post = fs.readFileSync(new URL("../src/post-analysis.tsx", import.meta.url), "utf8");
  const runtime = fs.readFileSync(new URL("../src/alivibe-msa-runtime.ts", import.meta.url), "utf8");
  const tree = fs.readFileSync(new URL("../src/lineage-tree-viewer.tsx", import.meta.url), "utf8");
  const styles = fs.readFileSync(new URL("../src/globals.css", import.meta.url), "utf8");
  assert.match(post, /MSA lineage · Alivibe WASM/);
  assert.match(post, /AA MSA · Alivibe WASM/);
  assert.match(post, /runDirectAlivibeAaMsa/);
  assert.match(post, /Drop non-productive/);
  assert.match(post, /workbenchLineageRows/);
  assert.match(post, /translateCodonSequence/);
  assert.match(post, /projectCodonAlignment/);
  assert.match(post, /runDirectAlivibeMsa/);
  assert.match(post, /runAlivibeMsaTask/);
  assert.match(post, /record\.sequence\.replaceAll\("-", ""\)/);
  assert.match(post, /Open \+ load in Alivibe/);
  assert.match(post, /Read live Alivibe NT view/);
  assert.match(post, /searchParams\.set\("swigBridge", String\(ALIVIBE_BRIDGE_VERSION\)\)/);
  assert.match(post, /searchParams\.set\("source", ALIVIBE_SOURCE_REVISION\.slice\(0, 12\)\)/);
  assert.match(runtime, /signal\?\.addEventListener\("abort", abort/);
  assert.match(runtime, /job\.cancel\(\)/);
  assert.match(runtime, /scoringMode: AlivibeMsaScoringMode = "nucleotide"/);
  assert.match(post, /Show sequence names/);
  assert.match(tree, /Show tip names/);
  assert.match(tree, /const labelWidth = showNames \? 208 : 12/);
  assert.match(styles, /\.lineage-alignment-preview\.names-hidden > div \{ grid-template-columns: max-content; \}/);
});

test("phylogenetic UCA posterior uses contour-bounded embedded glyphs and aligned HMM tracks", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const component = fs.readFileSync(new URL("../src/probability-logo.tsx", import.meta.url), "utf8");
  const glyphs = fs.readFileSync(new URL("../src/logo-glyphs.ts", import.meta.url), "utf8");
  const panel = fs.readFileSync(new URL("../src/phylo-uca/panel.tsx", import.meta.url), "utf8");
  const annotation = fs.readFileSync(new URL("../src/phylo-uca/hmm-annotation.tsx", import.meta.url), "utf8");
  const annotationModel = fs.readFileSync(new URL("../src/phylo-uca/hmm-annotation-model.ts", import.meta.url), "utf8");
  const hmm = fs.readFileSync(new URL("../src/phylo-uca/hmm.ts", import.meta.url), "utf8");
  const styles = fs.readFileSync(new URL("../src/globals.css", import.meta.url), "utf8");
  assert.match(component, /LOGO_MONOSPACE_FONT/);
  assert.match(component, /logoGlyphRun/);
  assert.match(component, /<path/);
  assert.match(component, /run\.yMax \* scaleY/);
  assert.match(component, /\$\{-scaleY\}/);
  assert.doesNotMatch(component, /getBBox|useLogoGlyphBounds|data-logo-measurement/);
  assert.match(glyphs, /Literal glyph contours from DejaVu Sans Mono Bold/);
  assert.match(glyphs, /"Q": \{[^\n]+yMin: -281/);
  assert.match(component, /width=\{columnWidth\}/);
  assert.doesNotMatch(component, /textLength=/);
  assert.match(component, /letter height is marginal frequency and is not scaled by entropy/i);
  assert.match(panel, /UCA posterior frequency logo/);
  assert.match(panel, /Logo SVG ↓/);
  assert.match(panel, /Full tracks \+ logo SVG ↓/);
  assert.match(panel, /Visible tracks \+ logo SVG ↓/);
  assert.match(panel, /serializePhyloUcaTrackLogoSvg/);
  assert.match(panel, /bottomAnnotations=\{logoBottomAnnotations\}/);
  assert.match(panel, /PhyloUcaPlacementMap/);
  assert.match(panel, /alignmentRetainedRows/);
  assert.match(panel, /retainedLineageRows\.map/);
  assert.match(panel, /V\/J nucleotide mixture · default/);
  assert.match(panel, /"Full-HMM edges"[\s\S]{0,120}?\(0 = all\)/);
  assert.match(panel, /Conditional ML · fastest/);
  assert.match(panel, /Explicit grid marginalization/);
  assert.match(panel, /Gibbs\/MH · default · continuous placement/);
  assert.match(panel, /Exact pendant-length grid/);
  assert.match(panel, /does not use the grid settings/);
  assert.match(panel, /Gibbs\/MH mixing/);
  assert.match(panel, />Best path</);
  assert.match(panel, />Marginalized</);
  assert.match(panel, /PhyloUcaHmmAnnotationTracks/);
  assert.match(panel, /collapseAndOrderHmmAnnotationTracks/);
  assert.match(panel, /labelOffset=\{annotationScrollLeft\}/);
  assert.match(panel, /Reset all UCA settings to defaults/);
  assert.match(panel, /Additional-D probability[\s\S]{0,500}?max="1"/);
  assert.doesNotMatch(panel, /Additional-D probability[\s\S]{0,500}?max="0\.5"/);
  assert.match(panel, /adoptedResultSnapshotRef\.current === incomingGeneratedAt/);
  assert.match(panel, /adoptedResultSnapshotRef\.current = inference\.generatedAt/);
  assert.doesNotMatch(panel, /\[matchingInitial, result\?\.generatedAt\]/);
  assert.match(app, /addEventListener\("beforeunload", warnBeforeLeaving\)/);
  assert.match(app, /addEventListener\("popstate", interceptHistoryDeparture\)/);
  assert.match(app, /setLeavePrompt\(true\)/);
  assert.match(app, /Stay on this page/);
  assert.match(app, /Leave anyway/);
  assert.match(styles, /overscroll-behavior-x:\s*none/);
  assert.match(panel, /position \* logoColumnWidth \/ 3/);
  assert.match(panel, />Codon</);
  assert.match(panel, /Amino acid/);
  assert.match(panel, /does not multiply nucleotide marginals/);
  assert.match(annotation, /total glyph height is the posterior track occupancy/i);
  assert.match(annotation, /FittedLogoGlyph/);
  assert.match(annotation, /nucleotide mass/);
  assert.match(annotation, /phylo-uca-track-label-panel/);
  assert.match(annotation, /exact scrolled crop/i);
  assert.match(annotation, /viewBox.*\$\{x\}.*\$\{y\}.*\$\{width\}.*\$\{height\}/);
  assert.match(annotation, /root\.append\(logo\)/);
  assert.match(annotation, /trackHeight \+ logoSize\.height/);
  assert.match(annotationModel, /combinedTrack\(`display\|\$\{key\}`/);
  assert.match(annotationModel, /weightedCenter/);
  assert.match(annotationModel, /mode === "marginalized"/);
  assert.match(annotationModel, /NT · all non-template mass/);
  assert.match(styles, /\.phylo-uca-track-legend \.v/);
  assert.match(styles, /\.phylo-uca-track-label-panel/);
  assert.match(hmm, /alignment site minus D-reference position/);
  assert.match(hmm, /Exact P\(x_i,x_\{i\+1\},x_\{i\+2\}/);

  for (const symbol of new Set([...NUCLEOTIDE_LOGO_SYMBOLS, ...AMINO_ACID_LOGO_SYMBOLS])) {
    const run = logoGlyphRun(symbol);
    assert.ok(run, `missing embedded outline for ${symbol}`);
    assert.equal(run.width, run.xMax - run.xMin);
    assert.equal(run.height, run.yMax - run.yMin);
    assert.ok(run.width > 0 && run.height > 0, `degenerate embedded outline for ${symbol}`);
  }
  const codon = logoGlyphRun("A-Q");
  assert.ok(codon && codon.paths.length === 3, "multi-character codon glyph runs must retain all contours");
});

test("a run is refused before it starts when an option lacks the germlines it needs", () => {
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const requestRun = app.slice(app.indexOf("function requestRun()"), app.indexOf("if (session && !window.confirm(", app.indexOf("function requestRun()")));
  assert.ok(requestRun.length > 0);
  // The engine would refuse both with its command-line wording; Swig Web says it first in its own.
  assert.match(requestRun, /if \(constantPrefixIdentityEnabled && !compiled\.counts\.C\) \{\s*setRunError\("The covered C-prefix identity gate requires at least one C germline in the composed reference set\."\);\s*return;/);
  assert.match(requestRun, /if \(!compiled\.counts\.D\) \{\s*setRunError\("Double-D screening requires at least one D germline in the composed reference set\."\);\s*return;/);
});

test("Swig Web shows the aligner's warnings in the run summary", () => {
  // Regression: the scheduler dropped the warnings the compute workers' engines return.
  const app = fs.readFileSync(new URL("../src/swig-app.tsx", import.meta.url), "utf8");
  const scheduler = fs.readFileSync(new URL("../src/aligner-worker.ts", import.meta.url), "utf8");
  const runtime = fs.readFileSync(new URL("../src/aligner-runtime.ts", import.meta.url), "utf8");
  assert.match(scheduler, /slot\.warnings = message\.warnings \?\? \[\]/);
  assert.match(scheduler, /warnings: slots\[0\]\?\.warnings \?\? \[\]/);
  assert.match(runtime, /warnings: message\.warnings \?\? \[\]/);
  assert.match(app, /for \(const warning of completed\.warnings\) alignerWarnings\.add\(warning\)/);
  assert.match(app, /alignerWarnings: result\.alignerWarnings/);
  // Both result layouts list them.
  assert.equal(app.match(/\{alignerNotes\}/g)?.length, 2);
  assert.match(app, /className="aligner-warnings" role="status"/);
});
