//#region src/aligner-engine.ts
const CALLING_PROFILES = [
	"r_optimized",
	"truth_optimized",
	"igblast_compatible",
	"sensitive_d"
];
function isCallingProfile(value) {
	return typeof value === "string" && CALLING_PROFILES.includes(value);
}
/** The settings as the engine's options, by the names of the native `swig` command line. */
function alignerOptions(settings) {
	const options = [
		["profile", settings.callingProfile],
		["strand", settings.strand === 1 ? "forward" : settings.strand === 2 ? "reverse" : "both"],
		["min-identity", String(settings.minimumIdentity)]
	];
	if (settings.minimumConstantPrefixIdentity != null) options.push(["c-prefix-identity", String(settings.minimumConstantPrefixIdentity)]);
	const doubleD = settings.doubleD;
	if (doubleD && doubleD.mode !== "off") options.push(["double-d", doubleD.mode], ["double-d-min-span", String(Math.round(doubleD.minimumVjSpan))], ["double-d-seed-length", String(Math.round(doubleD.seedLength))], ["double-d-pseudo-trim", String(Math.round(doubleD.pseudoTrim))], ["double-d-pseudo-mismatches", String(Math.round(doubleD.maximumPseudoMismatches))], ["double-d-min-gain", String(Math.round(doubleD.minimumScoreGain))]);
	const tuning = settings.tuning ?? {};
	if (tuning.minimumDMatch != null) options.push(["min-d-match", String(Math.round(tuning.minimumDMatch))]);
	if (tuning.dMismatch != null) options.push(["d-mismatch", String(Math.round(tuning.dMismatch))]);
	if (tuning.jMismatch != null) options.push(["j-mismatch", String(Math.round(tuning.jMismatch))]);
	if (tuning.minimumJLength != null) options.push(["min-j-length", String(Math.round(tuning.minimumJLength))]);
	return options;
}
const encoder = new TextEncoder();
const decoder = new TextDecoder();
/** Split `table` (header line, then rows) into its header and rows. */
function splitTable(table, what) {
	const newline = table.indexOf(10);
	if (newline < 0) throw new Error(`The Swig aligner returned an invalid ${what}.`);
	return {
		header: decoder.decode(table.subarray(0, newline)).replace(/\r$/, ""),
		body: table.slice(newline + 1)
	};
}
var AlignerEngine = class AlignerEngine {
	/** The database summary and warnings from loading the germlines. */
	summary;
	exports;
	constructor(exports, summary) {
		this.exports = exports;
		this.summary = summary;
	}
	/** One engine over `germlines` with `settings`, from the compiled module. */
	static async create(module, settings, germlines) {
		const exports = (await WebAssembly.instantiate(module, {})).exports;
		const call = (action) => {
			try {
				return action();
			} catch (error) {
				exports.swig_panic_message();
				const message = AlignerEngine.read(exports, exports.swig_error_ptr(), exports.swig_error_len());
				throw new Error(message || (error instanceof Error ? error.message : String(error)));
			}
		};
		for (const [name, value] of alignerOptions(settings)) {
			const nameBytes = AlignerEngine.put(exports, encoder.encode(name));
			const valueBytes = AlignerEngine.put(exports, encoder.encode(value));
			try {
				if (call(() => exports.swig_set_option(...nameBytes, ...valueBytes)) !== 0) throw new Error(AlignerEngine.error(exports));
			} finally {
				exports.swig_free(...nameBytes);
				exports.swig_free(...valueBytes);
			}
		}
		const texts = [
			germlines.V,
			germlines.D,
			germlines.J,
			germlines.C
		].map((text) => AlignerEngine.put(exports, encoder.encode(text || "")));
		try {
			const records = call(() => exports.swig_init_database(...texts.flat()));
			if (records < 0) throw new Error(AlignerEngine.error(exports));
			const summary = JSON.parse(AlignerEngine.read(exports, exports.swig_result_ptr(), exports.swig_result_len()));
			return {
				engine: new AlignerEngine(exports, summary),
				records
			};
		} finally {
			texts.forEach((buffer) => exports.swig_free(...buffer));
		}
	}
	/** The engine's version and alignment kernel. */
	version() {
		this.exports.swig_version();
		return JSON.parse(AlignerEngine.read(this.exports, this.exports.swig_result_ptr(), this.exports.swig_result_len()));
	}
	/**
	* Assign `count` reads (FASTA, FASTQ or an AIRR table with its header line):
	* the AIRR TSV header and one row per read in input order, and the double-D
	* sidecar when the screen is on.
	*/
	annotate(reads, count, doubleD) {
		const exports = this.exports;
		const input = AlignerEngine.put(exports, reads);
		let returned;
		try {
			returned = exports.swig_annotate(...input);
		} catch (error) {
			exports.swig_panic_message();
			throw new Error(AlignerEngine.error(exports) || (error instanceof Error ? error.message : String(error)));
		} finally {
			exports.swig_free(...input);
		}
		if (returned < 0) throw new Error(AlignerEngine.error(exports));
		if (returned !== count) throw new Error(`The Swig aligner returned ${returned} records for a ${count}-record input batch.`);
		const batch = {
			...splitTable(new Uint8Array(exports.memory.buffer, exports.swig_result_ptr(), exports.swig_result_len()), "AIRR table"),
			count
		};
		if (doubleD) {
			const { header, body } = splitTable(new Uint8Array(exports.memory.buffer, exports.swig_double_d_result_ptr(), exports.swig_double_d_result_len()), "double-D evidence table");
			batch.doubleDHeader = header;
			batch.doubleDBody = body;
			batch.doubleDCount = exports.swig_double_d_count();
		}
		return batch;
	}
	static put(exports, bytes) {
		const pointer = exports.swig_alloc(bytes.byteLength);
		if (!pointer) throw new Error("The Swig aligner ran out of WebAssembly memory.");
		new Uint8Array(exports.memory.buffer, pointer, bytes.byteLength).set(bytes);
		return [pointer, bytes.byteLength];
	}
	static read(exports, pointer, length) {
		return decoder.decode(new Uint8Array(exports.memory.buffer, pointer, length));
	}
	static error(exports) {
		return AlignerEngine.read(exports, exports.swig_error_ptr(), exports.swig_error_len()) || "The Swig aligner could not complete the request.";
	}
};
//#endregion
export { CALLING_PROFILES as n, isCallingProfile as r, AlignerEngine as t };
