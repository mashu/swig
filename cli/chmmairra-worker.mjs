import { d as prepareReferenceMsa, f as runChmm, g as threadSequenceToMsa, i as chmmairraDistanceFromReference, s as extendChmmAlignment } from "./chunks/post-analysis-core-B8-vIT7N.mjs";
import { parentPort } from "node:worker_threads";
//#region cli-src/chmmairra-worker.mjs
let msa = null;
let options = null;
let minDfr = 1;
const cache = /* @__PURE__ */ new Map();
const send = (message) => parentPort ? parentPort.postMessage(message) : globalThis.postMessage(message);
function initialize(message) {
	msa = prepareReferenceMsa(message.msa);
	options = message.options;
	minDfr = message.minDfr;
	cache.clear();
	return {
		references: msa.names.length,
		length: msa.length
	};
}
function batch(rows) {
	if (!msa || !options) throw new Error("CHMMAIRRa worker is not initialized.");
	const results = [];
	for (const input of rows) {
		if (!input.call || !input.sequenceAlignment || !input.germlineAlignment) {
			results.push({
				ordinal: input.ordinal,
				probability: NaN,
				dfr: 0,
				startingReference: "",
				recombinations: [],
				status: "missing_alignment"
			});
			continue;
		}
		const row = input.segment && input.flank ? {
			...input,
			...extendChmmAlignment(input.segment, input.sequenceAlignment, input.germlineAlignment, input.call, input.flank, msa)
		} : input;
		const dfr = chmmairraDistanceFromReference(row.sequenceAlignment, row.germlineAlignment);
		if (dfr < minDfr) {
			results.push({
				ordinal: row.ordinal,
				probability: NaN,
				dfr,
				startingReference: "",
				recombinations: [],
				status: "low_dfr"
			});
			continue;
		}
		const key = `${row.call}\0${row.sequenceAlignment}\0${row.germlineAlignment}`;
		const cached = cache.get(key);
		if (cached) {
			results.push({
				ordinal: row.ordinal,
				...cached
			});
			continue;
		}
		try {
			const threaded = threadSequenceToMsa(row.sequenceAlignment, row.germlineAlignment, row.call, msa);
			const value = {
				...runChmm(msa, threaded, row.sequenceAlignment, row.germlineAlignment, options),
				status: "evaluated"
			};
			cache.set(key, value);
			results.push({
				ordinal: row.ordinal,
				...value
			});
		} catch (error) {
			const value = {
				probability: NaN,
				dfr,
				startingReference: "",
				recombinations: [],
				status: "error",
				error: error instanceof Error ? error.message : String(error)
			};
			cache.set(key, value);
			results.push({
				ordinal: row.ordinal,
				...value
			});
		}
	}
	return { results };
}
const receive = (message) => {
	try {
		send({
			id: message.id,
			result: message.type === "init" ? initialize(message) : batch(message.rows)
		});
	} catch (error) {
		send({
			id: message.id,
			error: error instanceof Error ? error.message : String(error)
		});
	}
};
if (parentPort) parentPort.on("message", receive);
else globalThis.addEventListener("message", (event) => receive(event.data));
//#endregion
