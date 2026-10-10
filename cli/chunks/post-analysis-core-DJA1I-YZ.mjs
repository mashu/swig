//#region src/study-design.ts
const DEFAULT_PIPELINE_PLAN = {
	enabled: false,
	collapse: {
		enabled: true,
		mode: "exact",
		key: "trimmed",
		scope: "sample",
		unresolvedPolicy: "discard",
		respectConstantCall: true
	},
	chimera: {
		enabled: false,
		segment: "V",
		model: "auto",
		posteriorThreshold: .95,
		retainUnevaluated: true,
		msaSource: "selected",
		uploadedMsa: "",
		uploadedMsaName: ""
	},
	selection: {
		enabled: false,
		datasetId: "",
		sampleId: "",
		subjectId: "",
		cohort: "",
		timepoint: "",
		compartment: "",
		locus: "",
		vCall: "",
		vCallIncludeAmbiguous: false,
		jCall: "",
		jCallIncludeAmbiguous: false,
		cdr3Nt: "",
		cdr3Aa: "",
		productive: "any",
		hasCdr3: "any",
		doubleD: "any"
	},
	alleleRefinement: {
		enabled: false,
		model: "dirichlet",
		scope: "subject",
		segments: ["V", "J"],
		weighting: "unique",
		baselineNeighbourOdds: .01,
		shmLeakageSensitivity: 1,
		reassignmentPolicy: "confidence",
		applyMinimumPosterior: .8
	},
	lineage: {
		enabled: true,
		scope: "subject",
		identity: .85,
		resolution: "gene",
		ambiguity: "overlap",
		productiveOnly: true
	},
	shm: {
		enabled: true,
		metric: "vNtRate"
	},
	missingAlleles: { enabled: false }
};
function datasetScopeKey(record, scope = "global") {
	return `${scope}:${datasetScopeValue(record, scope)}`;
}
function datasetScopeValue(record, scope = "global") {
	if (scope === "global") return "complete study";
	if (scope === "dataset") return record.datasetId || "legacy";
	if (scope === "sample") return record.sampleId || record.datasetId || "legacy";
	if (scope === "subject") return record.subjectId || record.sampleId || record.datasetId || "legacy";
	return record.cohort || record.subjectId || record.sampleId || record.datasetId || "legacy";
}
const METADATA_FIELDS = [
	"swig_dataset_id",
	"sample_id",
	"subject_id",
	"swig_cohort",
	"swig_timepoint",
	"swig_compartment",
	"swig_source_sequence_id"
];
function normalizedLines(body) {
	return (typeof body === "string" ? body : new TextDecoder().decode(body)).split("\n").map((line) => line.replace(/\r$/, "")).filter(Boolean);
}
function replaceOrAppend(values, positions, field, value) {
	const position = positions.get(field);
	if (position === void 0) values.push(value);
	else values[position] = value;
}
/** Adds study metadata to an aligner AIRR batch without materializing the input file. */
function annotateAirrBatch(headerLine, body, dataset) {
	const headers = headerLine.replace(/\r$/, "").split("	");
	const positions = new Map(headers.map((field, index) => [field, index]));
	for (const field of METADATA_FIELDS) if (!positions.has(field)) headers.push(field);
	const sourceSequencePosition = positions.get("sequence_id");
	const bodyText = normalizedLines(body).map((line) => {
		const values = line.split("	");
		const sourceSequenceId = sourceSequencePosition === void 0 ? "" : values[sourceSequencePosition] ?? "";
		if (sourceSequencePosition !== void 0) values[sourceSequencePosition] = `${dataset.datasetId}::${sourceSequenceId || "record"}`;
		replaceOrAppend(values, positions, "swig_dataset_id", dataset.datasetId);
		replaceOrAppend(values, positions, "sample_id", dataset.sampleId);
		replaceOrAppend(values, positions, "subject_id", dataset.subjectId);
		replaceOrAppend(values, positions, "swig_cohort", dataset.cohort);
		replaceOrAppend(values, positions, "swig_timepoint", dataset.timepoint ?? "");
		replaceOrAppend(values, positions, "swig_compartment", dataset.compartment ?? "");
		replaceOrAppend(values, positions, "swig_source_sequence_id", sourceSequenceId);
		return values.join("	");
	}).join("\n");
	return {
		header: headers.join("	"),
		body: bodyText ? `${bodyText}\n` : ""
	};
}
/** Mirrors metadata and collision-safe IDs into the sparse double-D table. */
function annotateDoubleDBatch(headerLine, body, dataset) {
	const headers = headerLine.replace(/\r$/, "").split("	");
	const positions = new Map(headers.map((field, index) => [field, index]));
	for (const field of METADATA_FIELDS.slice(0, 6)) if (!positions.has(field)) headers.push(field);
	const sequencePosition = positions.get("sequence_id");
	const bodyText = normalizedLines(body).map((line) => {
		const values = line.split("	");
		if (sequencePosition !== void 0) values[sequencePosition] = `${dataset.datasetId}::${values[sequencePosition] || "record"}`;
		replaceOrAppend(values, positions, "swig_dataset_id", dataset.datasetId);
		replaceOrAppend(values, positions, "sample_id", dataset.sampleId);
		replaceOrAppend(values, positions, "subject_id", dataset.subjectId);
		replaceOrAppend(values, positions, "swig_cohort", dataset.cohort);
		replaceOrAppend(values, positions, "swig_timepoint", dataset.timepoint ?? "");
		replaceOrAppend(values, positions, "swig_compartment", dataset.compartment ?? "");
		return values.join("	");
	}).join("\n");
	return {
		header: headers.join("	"),
		body: bodyText ? `${bodyText}\n` : ""
	};
}
function stableDatasetSeed(baseSeed, datasetIndex) {
	return Math.trunc(baseSeed) + Math.imul(datasetIndex + 1, 1831565813) | 0;
}
//#endregion
//#region src/post-analysis-core.ts
const HASH_SEEDS = [
	2166136261,
	2654435769,
	2246822507,
	3266489909,
	668265263,
	374761393,
	3550635116,
	4251993797
];
function mix32(value) {
	value ^= value >>> 16;
	value = Math.imul(value, 2146121005);
	value ^= value >>> 15;
	value = Math.imul(value, 2221713035);
	value ^= value >>> 16;
	return value >>> 0;
}
function hashSequence(value, seed = HASH_SEEDS[0]) {
	let hash = seed >>> 0;
	for (let index = 0; index < value.length; index += 1) {
		hash ^= value.charCodeAt(index);
		hash = Math.imul(hash, 16777619);
	}
	return mix32(hash);
}
function sequenceFingerprint(value) {
	const normalized = normalizeNt(value);
	return `${normalized.length}:${[
		0,
		2,
		4,
		6
	].map((index) => hashSequence(normalized, HASH_SEEDS[index]).toString(36)).join(":")}`;
}
function normalizeNt(value) {
	return value.toUpperCase().replaceAll("U", "T").replace(/[^ACGTN]/g, "");
}
function normalizeCall(call, resolution) {
	const clean = call.trim().toUpperCase();
	return resolution === "gene" ? clean.replace(/\*.*$/, "") : clean;
}
function callSet(value, resolution, ambiguity) {
	const values = value.split(",").map((call) => normalizeCall(call, resolution)).filter(Boolean);
	return ambiguity === "top" ? values.slice(0, 1) : [...new Set(values)].sort();
}
function callsCompatible(left, right, resolution, ambiguity) {
	const a = callSet(left, resolution, ambiguity);
	const b = callSet(right, resolution, ambiguity);
	if (!a.length || !b.length) return false;
	if (ambiguity === "strict") return a.length === b.length && a.every((value, index) => value === b[index]);
	const values = new Set(a);
	return b.some((value) => values.has(value));
}
function hammingDistanceWithin(left, right, maximum, ambiguousN = true) {
	if (left.length !== right.length) return maximum + 1;
	let distance = 0;
	for (let index = 0; index < left.length; index += 1) if (left[index] !== right[index] && (!ambiguousN || left[index] !== "N" && right[index] !== "N")) {
		distance += 1;
		if (distance > maximum) return distance;
	}
	return distance;
}
function constantCallPartition(record, enabled) {
	if (!enabled) return "";
	return `\u0000C:${callSet(record.cCall ?? "", "gene", "top")[0] ?? "__C_UNASSIGNED__"}`;
}
function dedupKey(record, key, scope, respectConstantCall) {
	const prefix = `${datasetScopeKey(record, scope)}${constantCallPartition(record, respectConstantCall)}\u0000`;
	if (key === "sequence") return `${prefix}${record.sequenceFingerprint}`;
	if (key === "trimmed") return `${prefix}${record.trimmedFingerprint}`;
	if (key === "cdr3") return `${prefix}${record.locus}\u0000${record.cdr3Nt}`;
	return `${prefix}${record.locus}\u0000${record.vCall}\u0000${record.jCall}\u0000${record.cdr3Nt}`;
}
function largestCountGroups(counts, limit = 100) {
	const heap = [];
	const worse = (left, right) => counts[left] < counts[right] || counts[left] === counts[right] && left > right;
	const rise = (start) => {
		let index = start;
		while (index > 0) {
			const parent = Math.floor((index - 1) / 2);
			if (!worse(heap[index], heap[parent])) break;
			[heap[index], heap[parent]] = [heap[parent], heap[index]];
			index = parent;
		}
	};
	const sink = () => {
		let index = 0;
		while (true) {
			const left = index * 2 + 1;
			const right = left + 1;
			let next = index;
			if (left < heap.length && worse(heap[left], heap[next])) next = left;
			if (right < heap.length && worse(heap[right], heap[next])) next = right;
			if (next === index) break;
			[heap[index], heap[next]] = [heap[next], heap[index]];
			index = next;
		}
	};
	for (let ordinal = 0; ordinal < counts.length; ordinal += 1) {
		if (counts[ordinal] <= 1) continue;
		if (heap.length < limit) {
			heap.push(ordinal);
			rise(heap.length - 1);
		} else if (worse(heap[0], ordinal)) {
			heap[0] = ordinal;
			sink();
		}
	}
	return heap.sort((left, right) => counts[right] - counts[left] || left - right).map((ordinal) => ({
		ordinal,
		count: counts[ordinal]
	}));
}
function hasUsableDedupKey(record, key) {
	if (key === "sequence") return Number(record.sequenceFingerprint.split(":", 1)[0]) > 0;
	if (key === "trimmed") return Number(record.trimmedFingerprint.split(":", 1)[0]) > 0;
	if (key === "cdr3") return Boolean(record.locus && record.cdr3Nt);
	return Boolean(record.locus && record.vCall && record.jCall && record.cdr3Nt);
}
function deduplicate(records, key, unresolvedPolicy = "discard", scope = "global", respectConstantCall = true) {
	const representatives = new Int32Array(records.length);
	representatives.fill(-1);
	const counts = new Uint32Array(records.length);
	const seen = /* @__PURE__ */ new Map();
	let inputAbundance = 0;
	let unresolvedRecords = 0;
	let retainedUnresolved = 0;
	for (let index = 0; index < records.length; index += 1) {
		const weight = Math.max(1, Math.floor(records[index].inputCount ?? 1));
		inputAbundance += weight;
		if (!hasUsableDedupKey(records[index], key)) {
			unresolvedRecords += 1;
			if (unresolvedPolicy === "retain") {
				representatives[index] = index;
				counts[index] = weight;
				retainedUnresolved += 1;
			}
			continue;
		}
		const value = dedupKey(records[index], key, scope, respectConstantCall);
		const previous = seen.get(value);
		if (previous === void 0) {
			seen.set(value, index);
			representatives[index] = index;
			counts[index] = weight;
		} else {
			representatives[index] = previous;
			counts[previous] += weight;
		}
	}
	const largestGroups = largestCountGroups(counts);
	const uniqueRecords = seen.size + retainedUnresolved;
	return {
		mode: "exact",
		key,
		algorithm: "Exact key collapse",
		inputRecords: records.length,
		inputAbundance,
		uniqueRecords,
		collapsedRecords: records.length - uniqueRecords,
		representatives,
		counts,
		largestGroups,
		partitions: 1,
		candidateComparisons: 0,
		indelMergedVariants: 0,
		substitutionMergedVariants: 0,
		excludedAmbiguous: 0,
		unresolvedRecords,
		warnings: unresolvedRecords ? [`${unresolvedRecords.toLocaleString()} records without a usable ${key} key were ${unresolvedPolicy === "retain" ? "retained unchanged" : "discarded from the downstream representative set"}.`] : []
	};
}
function createExactDedupPlan(records, key, unresolvedPolicy = "discard", scope = "global", respectConstantCall = true, requestedShards = 1) {
	const shardCount = Math.max(1, Math.min(Math.floor(requestedShards) || 1, records.length || 1));
	const shardOrdinals = Array.from({ length: shardCount }, () => []);
	const shardWeights = Array.from({ length: shardCount }, () => []);
	const shardKeys = Array.from({ length: shardCount }, () => []);
	const unresolvedOrdinals = [];
	const unresolvedWeights = [];
	let inputAbundance = 0;
	for (let ordinal = 0; ordinal < records.length; ordinal += 1) {
		const weight = Math.max(1, Math.floor(records[ordinal].inputCount ?? 1));
		inputAbundance += weight;
		if (!hasUsableDedupKey(records[ordinal], key)) {
			unresolvedOrdinals.push(ordinal);
			unresolvedWeights.push(weight);
			continue;
		}
		const value = dedupKey(records[ordinal], key, scope, respectConstantCall);
		const shard = hashSequence(value) % shardCount;
		shardOrdinals[shard].push(ordinal);
		shardWeights[shard].push(weight);
		shardKeys[shard].push(value);
	}
	const jobs = [];
	for (let shard = 0; shard < shardCount; shard += 1) {
		if (!shardOrdinals[shard].length) continue;
		jobs.push({
			id: shard,
			ordinals: Int32Array.from(shardOrdinals[shard]),
			weights: Uint32Array.from(shardWeights[shard]),
			keys: shardKeys[shard]
		});
	}
	return {
		key,
		recordCount: records.length,
		inputAbundance,
		unresolvedPolicy,
		unresolvedOrdinals: Int32Array.from(unresolvedOrdinals),
		unresolvedWeights: Uint32Array.from(unresolvedWeights),
		jobs
	};
}
function runExactDedupJob(job) {
	if (job.ordinals.length !== job.weights.length || job.ordinals.length !== job.keys.length) throw new Error("The exact-collapse worker received inconsistent shard vectors.");
	const representatives = new Int32Array(job.ordinals.length);
	const seen = /* @__PURE__ */ new Map();
	const countByRepresentative = /* @__PURE__ */ new Map();
	for (let index = 0; index < job.ordinals.length; index += 1) {
		const previous = seen.get(job.keys[index]);
		const representative = previous ?? job.ordinals[index];
		if (previous === void 0) seen.set(job.keys[index], representative);
		representatives[index] = representative;
		countByRepresentative.set(representative, (countByRepresentative.get(representative) ?? 0) + job.weights[index]);
	}
	const representativeOrdinals = Int32Array.from(countByRepresentative.keys());
	const representativeCounts = Uint32Array.from(representativeOrdinals, (ordinal) => countByRepresentative.get(ordinal) ?? 0);
	return {
		id: job.id,
		representatives,
		representativeOrdinals,
		representativeCounts,
		uniqueRecords: seen.size
	};
}
function finishExactDedupPlan(plan, results) {
	const representatives = new Int32Array(plan.recordCount);
	representatives.fill(-1);
	const counts = new Uint32Array(plan.recordCount);
	const resultById = new Map(results.map((result) => [result.id, result]));
	let uniqueRecords = 0;
	for (const job of plan.jobs) {
		const result = resultById.get(job.id);
		if (!result || result.representatives.length !== job.ordinals.length) throw new Error(`Exact-collapse shard ${job.id} did not return a complete result.`);
		uniqueRecords += result.uniqueRecords;
		for (let index = 0; index < job.ordinals.length; index += 1) representatives[job.ordinals[index]] = result.representatives[index];
		for (let index = 0; index < result.representativeOrdinals.length; index += 1) counts[result.representativeOrdinals[index]] = result.representativeCounts[index];
	}
	if (plan.unresolvedPolicy === "retain") {
		uniqueRecords += plan.unresolvedOrdinals.length;
		for (let index = 0; index < plan.unresolvedOrdinals.length; index += 1) {
			const ordinal = plan.unresolvedOrdinals[index];
			representatives[ordinal] = ordinal;
			counts[ordinal] = plan.unresolvedWeights[index];
		}
	}
	const unresolvedRecords = plan.unresolvedOrdinals.length;
	return {
		mode: "exact",
		key: plan.key,
		algorithm: "Exact key collapse",
		inputRecords: plan.recordCount,
		inputAbundance: plan.inputAbundance,
		uniqueRecords,
		collapsedRecords: plan.recordCount - uniqueRecords,
		representatives,
		counts,
		largestGroups: largestCountGroups(counts),
		partitions: 1,
		candidateComparisons: 0,
		indelMergedVariants: 0,
		substitutionMergedVariants: 0,
		excludedAmbiguous: 0,
		unresolvedRecords,
		warnings: unresolvedRecords ? [`${unresolvedRecords.toLocaleString()} records without a usable ${plan.key} key were ${plan.unresolvedPolicy === "retain" ? "retained unchanged" : "discarded from the downstream representative set"}.`] : []
	};
}
/**
* Append-only two-bit DNA storage. Denoising may need the complete trimmed VDJ
* sequence, but retaining one JavaScript string per unique read is prohibitive
* for large repertoires. Fixed chunks also avoid a transient 2x allocation
* when an expanding typed array is copied.
*/
var PackedDnaArena = class {
	chunkWords = 262144;
	chunks = [];
	used = [];
	append(sequence) {
		const words = Math.ceil(sequence.length / 16);
		let chunk = this.chunks.length - 1;
		if (chunk < 0 || this.used[chunk] + words > this.chunkWords) {
			chunk = this.chunks.length;
			this.chunks.push(new Uint32Array(Math.max(this.chunkWords, words)));
			this.used.push(0);
		}
		const offset = this.used[chunk];
		this.used[chunk] += words;
		const target = this.chunks[chunk];
		for (let index = 0; index < sequence.length; index += 1) {
			const code = sequence[index] === "C" ? 1 : sequence[index] === "G" ? 2 : sequence[index] === "T" ? 3 : 0;
			target[offset + (index >>> 4)] |= code << (index & 15) * 2;
		}
		return {
			chunk,
			offset,
			length: sequence.length
		};
	}
	base(location, index) {
		const code = this.chunks[location.chunk][location.offset + (index >>> 4)] >>> (index & 15) * 2 & 3;
		return code === 0 ? "A" : code === 1 ? "C" : code === 2 ? "G" : "T";
	}
	equals(location, sequence) {
		if (location.length !== sequence.length) return false;
		for (let index = 0; index < sequence.length; index += 1) if (this.base(location, index) !== sequence[index]) return false;
		return true;
	}
	decode(location) {
		let sequence = "";
		for (let start = 0; start < location.length; start += 1024) {
			const values = [];
			const end = Math.min(location.length, start + 1024);
			for (let index = start; index < end; index += 1) values.push(this.base(location, index));
			sequence += values.join("");
		}
		return sequence;
	}
};
function denoisePartition(record, options) {
	const v = callSet(record.vCall, options.callResolution, options.ambiguity);
	const j = callSet(record.jCall, options.callResolution, options.ambiguity);
	if (!v.length || !j.length) return null;
	return `${datasetScopeKey(record, options.scope ?? "global")}${constantCallPartition(record, options.respectConstantCall ?? true)}\u0000${record.locus || "?"}\u0000${v.join("+")}\u0000${j.join("+")}`;
}
function kmerProfile(sequence, blockCount, k = 6) {
	const space = 1 << 2 * k;
	const counts = new Uint16Array(space);
	const touched = [];
	const mask = space - 1;
	let code = 0;
	for (let index = 0; index < sequence.length; index += 1) {
		const base = sequence[index] === "C" ? 1 : sequence[index] === "G" ? 2 : sequence[index] === "T" ? 3 : 0;
		code = (code << 2 | base) & mask;
		if (index < k - 1) continue;
		if (!counts[code]) touched.push(code);
		if (counts[code] < 65535) counts[code] += 1;
	}
	touched.sort((left, right) => left - right);
	const codes = Uint16Array.from(touched);
	const sparseCounts = Uint16Array.from(touched, (value) => counts[value]);
	const hashes = new Uint32Array(blockCount);
	for (let block = 0; block < blockCount; block += 1) hashes[block] = mix32(2654435769 ^ block);
	for (let index = 0; index < codes.length; index += 1) {
		const block = codes[index] % blockCount;
		hashes[block] = mix32(hashes[block] ^ mix32(Math.imul(codes[index] + 1, 2246822507) ^ sparseCounts[index]));
	}
	return {
		codes,
		counts: sparseCounts,
		hashes
	};
}
function kmerSquaredDistance(left, right, ceiling = Number.POSITIVE_INFINITY) {
	let a = 0;
	let b = 0;
	let distance = 0;
	while (a < left.codes.length || b < right.codes.length) {
		if (b >= right.codes.length || a < left.codes.length && left.codes[a] < right.codes[b]) {
			distance += left.counts[a] ** 2;
			a += 1;
		} else if (a >= left.codes.length || right.codes[b] < left.codes[a]) {
			distance += right.counts[b] ** 2;
			b += 1;
		} else {
			distance += (left.counts[a] - right.counts[b]) ** 2;
			a += 1;
			b += 1;
		}
		if (distance > ceiling) return distance;
	}
	return distance;
}
/** Exact metric index for FAD's final nearest-centroid assignment. */
function buildKmerVpTree(points, profiles, compared) {
	if (!points.length) return null;
	const point = points[points.length - 1];
	if (points.length === 1) return {
		point,
		radius: 0,
		inside: null,
		outside: null
	};
	const distances = points.slice(0, -1).map((candidate) => {
		compared();
		return {
			candidate,
			distance: Math.sqrt(kmerSquaredDistance(profiles.get(point), profiles.get(candidate)))
		};
	}).sort((left, right) => left.distance - right.distance || left.candidate - right.candidate);
	const middle = Math.floor(distances.length / 2);
	return {
		point,
		radius: distances[middle]?.distance ?? 0,
		inside: buildKmerVpTree(distances.slice(0, middle).map((value) => value.candidate), profiles, compared),
		outside: buildKmerVpTree(distances.slice(middle).map((value) => value.candidate), profiles, compared)
	};
}
function nearestKmerPoint(tree, query, profiles, abundance, compared) {
	let best = tree.point;
	let bestDistance = Number.POSITIVE_INFINITY;
	const visit = (node) => {
		if (!node) return;
		compared();
		const distance = Math.sqrt(kmerSquaredDistance(profiles.get(query), profiles.get(node.point)));
		if (distance < bestDistance || distance === bestDistance && abundance(node.point) > abundance(best)) {
			best = node.point;
			bestDistance = distance;
		}
		if (!node.inside && !node.outside) return;
		if (distance < node.radius) {
			visit(node.inside);
			if (distance + bestDistance >= node.radius) visit(node.outside);
		} else {
			visit(node.outside);
			if (distance - bestDistance <= node.radius) visit(node.inside);
		}
	};
	visit(tree);
	return best;
}
function logGamma(value) {
	const coefficients = [
		.9999999999998099,
		676.5203681218851,
		-1259.1392167224028,
		771.3234287776531,
		-176.6150291621406,
		12.50734327868691,
		-.1385710952657201,
		9984369578019572e-21,
		1.505632735149312e-7
	];
	if (value < .5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * value)) - logGamma(1 - value);
	let x = coefficients[0];
	const z = value - 1;
	for (let index = 1; index < coefficients.length; index += 1) x += coefficients[index] / (z + index);
	const t = z + 7.5;
	return .5 * Math.log(2 * Math.PI) + (z + .5) * Math.log(t) - t + Math.log(x);
}
function regularizedGammaP(shape, value) {
	if (!(shape > 0) || value < 0 || Number.isNaN(value)) return NaN;
	if (value === 0) return 0;
	const logScale = -value + shape * Math.log(value) - logGamma(shape);
	if (value < shape + 1) {
		let term = 1 / shape;
		let sum = term;
		let denominator = shape;
		for (let iteration = 0; iteration < 1e4; iteration += 1) {
			denominator += 1;
			term *= value / denominator;
			sum += term;
			if (Math.abs(term) < Math.abs(sum) * 1e-14) break;
		}
		return Math.min(1, Math.max(0, sum * Math.exp(logScale)));
	}
	let b = value + 1 - shape;
	const floor = 1e-300;
	let c = 1 / floor;
	let d = 1 / Math.max(floor, b);
	let fraction = d;
	for (let iteration = 1; iteration < 1e4; iteration += 1) {
		const coefficient = -iteration * (iteration - shape);
		b += 2;
		d = coefficient * d + b;
		if (Math.abs(d) < floor) d = floor;
		c = b + coefficient / c;
		if (Math.abs(c) < floor) c = floor;
		d = 1 / d;
		const delta = d * c;
		fraction *= delta;
		if (Math.abs(delta - 1) < 1e-14) break;
	}
	const q = Math.exp(logScale) * fraction;
	return Math.min(1, Math.max(0, 1 - q));
}
/** P(Poisson(lambda) > observed), matching Distributions.jl ccdf. */
function poissonStrictUpperTail(observed, lambda) {
	if (!(lambda > 0)) return 0;
	if (observed < 0) return 1;
	return regularizedGammaP(Math.floor(observed) + 1, lambda);
}
function alternativeCount(length, distance) {
	if (distance <= 0) return 1;
	let combinations = 1;
	for (let index = 1; index <= distance; index += 1) combinations *= (length - distance + index) / index;
	return Math.max(1, combinations * 3 ** distance);
}
function sequenceBlocks(sequence, count) {
	const result = [];
	for (let index = 0; index < count; index += 1) {
		const start = Math.floor(index * sequence.length / count);
		const end = Math.floor((index + 1) * sequence.length / count);
		result.push(sequence.slice(start, end));
	}
	return result;
}
function indexedEditSegments(sequence, count) {
	const result = [];
	for (let index = 0; index < count; index += 1) {
		const start = Math.floor(index * sequence.length / count);
		const end = Math.floor((index + 1) * sequence.length / count);
		result.push({
			index,
			start,
			length: end - start,
			value: sequence.slice(start, end)
		});
	}
	return result;
}
/**
* Allocation-bounded Ukkonen-style edit profiler. The narrow diagonal band is
* O(length × maximum), and tie-breaking chooses fewer indels when two optimal
* paths exist so the aggressive indel rule is not triggered by an arbitrary
* alignment of equal-length sequences.
*/
function createBoundedEditProfiler(maximum) {
	const width = maximum * 2 + 3;
	const infinity = maximum + 1;
	let previousCost = new Int16Array(width);
	let previousSubstitutions = new Int16Array(width);
	let previousInsertions = new Int16Array(width);
	let previousDeletions = new Int16Array(width);
	let currentCost = new Int16Array(width);
	let currentSubstitutions = new Int16Array(width);
	let currentInsertions = new Int16Array(width);
	let currentDeletions = new Int16Array(width);
	return (parent, child) => {
		if (Math.abs(parent.length - child.length) > maximum) return null;
		previousCost.fill(infinity);
		previousSubstitutions.fill(0);
		previousInsertions.fill(0);
		previousDeletions.fill(0);
		for (let column = 0; column <= Math.min(child.length, maximum); column += 1) {
			const offset = column + maximum + 1;
			previousCost[offset] = column;
			previousInsertions[offset] = column;
		}
		for (let row = 1; row <= parent.length; row += 1) {
			currentCost.fill(infinity);
			currentSubstitutions.fill(0);
			currentInsertions.fill(0);
			currentDeletions.fill(0);
			const firstColumn = Math.max(0, row - maximum);
			const lastColumn = Math.min(child.length, row + maximum);
			for (let column = firstColumn; column <= lastColumn; column += 1) {
				const offset = column - row + maximum + 1;
				let bestCost = infinity;
				let bestSubstitutions = 0;
				let bestInsertions = 0;
				let bestDeletions = 0;
				const consider = (cost, substitutions, insertions, deletions) => {
					const indels = insertions + deletions;
					const bestIndels = bestInsertions + bestDeletions;
					if (cost < bestCost || cost === bestCost && indels < bestIndels || cost === bestCost && indels === bestIndels && substitutions < bestSubstitutions || cost === bestCost && indels === bestIndels && substitutions === bestSubstitutions && deletions < bestDeletions) {
						bestCost = cost;
						bestSubstitutions = substitutions;
						bestInsertions = insertions;
						bestDeletions = deletions;
					}
				};
				if (column > 0 && previousCost[offset] <= maximum) {
					const mismatch = parent[row - 1] === child[column - 1] ? 0 : 1;
					consider(previousCost[offset] + mismatch, previousSubstitutions[offset] + mismatch, previousInsertions[offset], previousDeletions[offset]);
				}
				if (previousCost[offset + 1] <= maximum) consider(previousCost[offset + 1] + 1, previousSubstitutions[offset + 1], previousInsertions[offset + 1], previousDeletions[offset + 1] + 1);
				if (column > 0 && currentCost[offset - 1] <= maximum) consider(currentCost[offset - 1] + 1, currentSubstitutions[offset - 1], currentInsertions[offset - 1] + 1, currentDeletions[offset - 1]);
				currentCost[offset] = bestCost;
				currentSubstitutions[offset] = bestSubstitutions;
				currentInsertions[offset] = bestInsertions;
				currentDeletions[offset] = bestDeletions;
			}
			[previousCost, currentCost] = [currentCost, previousCost];
			[previousSubstitutions, currentSubstitutions] = [currentSubstitutions, previousSubstitutions];
			[previousInsertions, currentInsertions] = [currentInsertions, previousInsertions];
			[previousDeletions, currentDeletions] = [currentDeletions, previousDeletions];
		}
		const offset = child.length - parent.length + maximum + 1;
		const distance = previousCost[offset];
		if (distance > maximum) return null;
		return {
			distance,
			substitutions: previousSubstitutions[offset],
			insertions: previousInsertions[offset],
			deletions: previousDeletions[offset]
		};
	};
}
/**
* Pure partition kernel shared by browser workers, CLI worker threads, and the
* serial fallback. A V/J partition is an exact independence boundary for all
* three denoisers; splitting inside it would change parent/centroid choices.
*/
function runDenoisePartitionJob(job) {
	const variants = job.variants;
	if (job.variantIndices.length !== variants.length) throw new Error(`Denoising partition ${job.id} has inconsistent variant vectors.`);
	const targets = Int32Array.from(variants, (_, index) => index);
	const ordered = Array.from({ length: variants.length }, (_, index) => index).sort((left, right) => variants[right].count - variants[left].count || variants[left].representative - variants[right].representative);
	const options = job.options;
	if (options.mode === "fad") {
		const maximumSquared = Math.max(0, Math.floor(12 * options.fadNeighborThreshold + 1e-9));
		const blockCount = Math.max(1, maximumSquared + 1);
		const profiles = /* @__PURE__ */ new Map();
		for (let index = 0; index < variants.length; index += 1) profiles.set(index, kmerProfile(variants[index].sequence, blockCount));
		const accepted = [];
		const blockIndex = /* @__PURE__ */ new Map();
		const candidates = /* @__PURE__ */ new Set();
		let comparisons = 0;
		let truncated = 0;
		const addAccepted = (variantIndex) => {
			accepted.push(variantIndex);
			profiles.get(variantIndex).hashes.forEach((hash, block) => {
				const key = `${block}:${hash}`;
				const values = blockIndex.get(key);
				if (values) values.push(variantIndex);
				else blockIndex.set(key, [variantIndex]);
			});
		};
		for (const variantIndex of ordered.filter((value) => variants[value].count >= options.minimumParentCount)) {
			const profile = profiles.get(variantIndex);
			candidates.clear();
			profile.hashes.forEach((hash, block) => {
				if (candidates.size >= options.maxCandidatesPerVariant) return;
				for (const candidate of blockIndex.get(`${block}:${hash}`) ?? []) {
					if (candidates.size >= options.maxCandidatesPerVariant) break;
					candidates.add(candidate);
				}
			});
			if (candidates.size >= options.maxCandidatesPerVariant) truncated += 1;
			const neighbors = [];
			for (const candidate of candidates) {
				comparisons += 1;
				const distance = kmerSquaredDistance(profile, profiles.get(candidate), maximumSquared);
				if (distance <= maximumSquared) neighbors.push({
					index: candidate,
					distance
				});
			}
			if (!neighbors.length) {
				addAccepted(variantIndex);
				continue;
			}
			neighbors.sort((left, right) => variants[right.index].count - variants[left.index].count || left.distance - right.distance || variants[left.index].representative - variants[right.index].representative);
			const parent = neighbors[0].index;
			const child = variants[variantIndex];
			const lambda = variants[parent].count / Math.max(Number.MIN_VALUE, options.expectedZeroErrorFraction) * options.errorRate;
			const adjusted = Math.min(1, poissonStrictUpperTail(child.count, lambda) * (child.sequence.length || 1));
			if (options.fadMethod === 2 && adjusted < options.alpha) addAccepted(variantIndex);
			else targets[variantIndex] = parent;
		}
		if (!accepted.length && ordered.length) addAccepted(ordered[0]);
		const acceptedSet = new Set(accepted);
		const vpTree = buildKmerVpTree(accepted, profiles, () => {
			comparisons += 1;
		});
		for (const variantIndex of ordered) targets[variantIndex] = acceptedSet.has(variantIndex) ? variantIndex : vpTree ? nearestKmerPoint(vpTree, variantIndex, profiles, (point) => variants[point].count, () => {
			comparisons += 1;
		}) : variantIndex;
		return {
			id: job.id,
			targets,
			comparisons,
			truncated,
			indelMergedVariants: 0,
			substitutionMergedVariants: 0,
			variantWork: variants.length + variants.reduce((count, variant) => count + (variant.count >= options.minimumParentCount ? 1 : 0), 0)
		};
	}
	if (options.mode === "indel") {
		const distanceLimit = options.maximumEditDistance;
		const blockCount = distanceLimit + 1;
		const parentIndex = /* @__PURE__ */ new Map();
		const shortParentsByLength = /* @__PURE__ */ new Map();
		const profileEdit = createBoundedEditProfiler(distanceLimit);
		let comparisons = 0;
		let truncated = 0;
		let indelMergedVariants = 0;
		let substitutionMergedVariants = 0;
		const candidates = /* @__PURE__ */ new Set();
		const addParent = (variantIndex) => {
			if (variants[variantIndex].count < options.minimumParentCount) return;
			const sequence = variants[variantIndex].sequence;
			const segments = indexedEditSegments(sequence, blockCount);
			if (segments.some((segment) => segment.length === 0)) {
				const values = shortParentsByLength.get(sequence.length);
				if (values) values.push(variantIndex);
				else shortParentsByLength.set(sequence.length, [variantIndex]);
				return;
			}
			for (const segment of segments) {
				const key = `${sequence.length}:${segment.index}:${segment.value}`;
				const values = parentIndex.get(key);
				if (values) values.push(variantIndex);
				else parentIndex.set(key, [variantIndex]);
			}
		};
		for (const variantIndex of ordered) {
			const child = variants[variantIndex];
			const sequence = child.sequence;
			candidates.clear();
			let capped = false;
			const addCandidates = (values) => {
				for (const candidate of values) {
					if (candidates.has(candidate)) continue;
					if (candidates.size >= options.maxCandidatesPerVariant) {
						capped = true;
						return;
					}
					candidates.add(candidate);
				}
			};
			const minimumParentLength = Math.max(1, sequence.length - distanceLimit);
			const maximumParentLength = sequence.length + distanceLimit;
			for (let parentLength = minimumParentLength; parentLength <= maximumParentLength && !capped; parentLength += 1) {
				addCandidates(shortParentsByLength.get(parentLength) ?? []);
				for (let segmentIndex = 0; segmentIndex < blockCount && !capped; segmentIndex += 1) {
					const parentStart = Math.floor(segmentIndex * parentLength / blockCount);
					const segmentLength = Math.floor((segmentIndex + 1) * parentLength / blockCount) - parentStart;
					if (!segmentLength) continue;
					const firstStart = Math.max(0, parentStart - distanceLimit);
					const lastStart = Math.min(sequence.length - segmentLength, parentStart + distanceLimit);
					for (let queryStart = firstStart; queryStart <= lastStart && !capped; queryStart += 1) addCandidates(parentIndex.get(`${parentLength}:${segmentIndex}:${sequence.slice(queryStart, queryStart + segmentLength)}`) ?? []);
				}
			}
			if (capped) truncated += 1;
			let best = -1;
			let bestProfile = null;
			for (const candidate of candidates) {
				const parent = variants[candidate];
				if (parent.count <= child.count) continue;
				comparisons += 1;
				const profile = profileEdit(parent.sequence, sequence);
				if (!profile || profile.distance < 1) continue;
				const indels = profile.insertions + profile.deletions;
				let plausible = false;
				if (indels > 0) plausible = parent.count / child.count >= options.minimumIndelParentRatio;
				else {
					const exactErrorProbability = (options.errorRate / 3) ** profile.substitutions * (1 - options.errorRate) ** Math.max(0, sequence.length - profile.substitutions);
					const lambda = parent.count * exactErrorProbability;
					plausible = Math.min(1, poissonStrictUpperTail(child.count, lambda) * alternativeCount(sequence.length, profile.substitutions)) >= options.alpha;
				}
				if (!plausible) continue;
				const bestParent = best >= 0 ? variants[best] : null;
				if (!bestProfile || profile.distance < bestProfile.distance || profile.distance === bestProfile.distance && profile.substitutions < bestProfile.substitutions || profile.distance === bestProfile.distance && profile.substitutions === bestProfile.substitutions && parent.count > (bestParent?.count ?? -1) || profile.distance === bestProfile.distance && profile.substitutions === bestProfile.substitutions && parent.count === bestParent?.count && parent.representative < (bestParent?.representative ?? Number.POSITIVE_INFINITY)) {
					best = candidate;
					bestProfile = profile;
				}
			}
			if (best >= 0 && bestProfile) {
				targets[variantIndex] = best;
				if (bestProfile.insertions + bestProfile.deletions > 0) indelMergedVariants += 1;
				else substitutionMergedVariants += 1;
			} else {
				targets[variantIndex] = variantIndex;
				addParent(variantIndex);
			}
		}
		return {
			id: job.id,
			targets,
			comparisons,
			truncated,
			indelMergedVariants,
			substitutionMergedVariants,
			variantWork: variants.length
		};
	}
	const distanceLimit = options.maximumHammingDistance;
	const blockCount = distanceLimit + 1;
	const parentIndex = /* @__PURE__ */ new Map();
	let comparisons = 0;
	let truncated = 0;
	let substitutionMergedVariants = 0;
	const candidates = /* @__PURE__ */ new Set();
	const addParent = (variantIndex) => {
		if (variants[variantIndex].count < options.minimumParentCount) return;
		const sequence = variants[variantIndex].sequence;
		sequenceBlocks(sequence, blockCount).forEach((block, blockIndex) => {
			const key = `${sequence.length}:${blockIndex}:${block}`;
			const values = parentIndex.get(key);
			if (values) values.push(variantIndex);
			else parentIndex.set(key, [variantIndex]);
		});
	};
	for (const variantIndex of ordered) {
		const child = variants[variantIndex];
		const sequence = child.sequence;
		candidates.clear();
		sequenceBlocks(sequence, blockCount).forEach((block, blockIndex) => {
			if (candidates.size >= options.maxCandidatesPerVariant) return;
			for (const candidate of parentIndex.get(`${sequence.length}:${blockIndex}:${block}`) ?? []) {
				if (candidates.size >= options.maxCandidatesPerVariant) break;
				candidates.add(candidate);
			}
		});
		if (candidates.size >= options.maxCandidatesPerVariant) truncated += 1;
		let best = -1;
		let bestLambda = -1;
		let bestDistance = Number.POSITIVE_INFINITY;
		for (const candidate of candidates) {
			comparisons += 1;
			const distance = hammingDistanceWithin(sequence, variants[candidate].sequence, distanceLimit, false);
			if (distance < 1 || distance > distanceLimit) continue;
			const lambda = variants[candidate].count * ((options.errorRate / 3) ** distance * (1 - options.errorRate) ** Math.max(0, sequence.length - distance));
			if (Math.min(1, poissonStrictUpperTail(child.count, lambda) * alternativeCount(sequence.length, distance)) >= options.alpha && (lambda > bestLambda || lambda === bestLambda && distance < bestDistance)) {
				best = candidate;
				bestLambda = lambda;
				bestDistance = distance;
			}
		}
		if (best >= 0) {
			targets[variantIndex] = best;
			substitutionMergedVariants += 1;
		} else {
			targets[variantIndex] = variantIndex;
			addParent(variantIndex);
		}
	}
	return {
		id: job.id,
		targets,
		comparisons,
		truncated,
		indelMergedVariants: 0,
		substitutionMergedVariants,
		variantWork: variants.length
	};
}
/**
* Streaming builder used by the post-analysis worker. Exact variants are
* dereplicated while batches are scanned, sequences live in a compact 2-bit
* arena, and only one V/J partition's temporary neighbor index exists at a
* time during finalization.
*/
var DenoiseAccumulator = class {
	arena = new PackedDnaArena();
	variants = [];
	variantsByFingerprint = /* @__PURE__ */ new Map();
	variantByOrdinal;
	processed;
	standalone = [];
	retainedAmbiguous = /* @__PURE__ */ new Map();
	excludedAmbiguous = 0;
	unresolvedRecords = 0;
	records;
	options;
	partitionJobsCache = null;
	constructor(records, options) {
		this.records = records;
		this.options = options;
		this.variantByOrdinal = new Int32Array(records.length);
		this.variantByOrdinal.fill(-1);
		this.processed = new Uint8Array(records.length);
		if (!(options.errorRate > 0 && options.errorRate < 1)) throw new Error("The denoising error rate must be between 0 and 1.");
		if (!(options.alpha > 0 && options.alpha < 1)) throw new Error("The denoising alpha must be between 0 and 1.");
		if (options.fadNeighborThreshold < 0) throw new Error("The FAD neighbor threshold cannot be negative.");
		if (options.maximumHammingDistance < 1 || options.maximumHammingDistance > 4) throw new Error("The conservative Hamming radius must be from 1 to 4.");
		if (!Number.isInteger(options.maximumEditDistance) || options.maximumEditDistance < 1 || options.maximumEditDistance > 2) throw new Error("The indel-aware edit radius must be 1 or 2.");
		if (!(options.minimumIndelParentRatio > 1)) throw new Error("The indel parent:child abundance ratio must be greater than 1.");
	}
	add(ordinal, rawSequence) {
		if (this.partitionJobsCache) throw new Error("Cannot add reads after denoising partition finalization has started.");
		if (ordinal < 0 || ordinal >= this.records.length || this.processed[ordinal]) return;
		this.processed[ordinal] = 1;
		const record = this.records[ordinal];
		const sequence = normalizeNt(rawSequence);
		const partition = denoisePartition(record, this.options);
		if (!sequence || !partition) {
			if (this.options.unresolvedPolicy === "retain") this.standalone.push(ordinal);
			this.unresolvedRecords += 1;
			return;
		}
		if (sequence.includes("N")) {
			if (this.options.ambiguousPolicy === "exclude") {
				this.excludedAmbiguous += 1;
				return;
			}
			const key = `${partition}\u0000${sequence}`;
			const existing = this.retainedAmbiguous.get(key);
			const weight = Math.max(1, Math.floor(record.inputCount ?? 1));
			if (existing) {
				existing.members.push(ordinal);
				existing.count += weight;
			} else this.retainedAmbiguous.set(key, {
				representative: ordinal,
				members: [ordinal],
				count: weight
			});
			return;
		}
		const fingerprint = `${partition}\u0000${sequenceFingerprint(sequence)}`;
		const candidates = this.variantsByFingerprint.get(fingerprint) ?? [];
		let variantIndex = candidates.find((candidate) => this.arena.equals(this.variants[candidate].location, sequence));
		const weight = Math.max(1, Math.floor(record.inputCount ?? 1));
		if (variantIndex === void 0) {
			variantIndex = this.variants.length;
			this.variants.push({
				location: this.arena.append(sequence),
				partition,
				representative: ordinal,
				count: weight,
				target: variantIndex
			});
			candidates.push(variantIndex);
			this.variantsByFingerprint.set(fingerprint, candidates);
		} else this.variants[variantIndex].count += weight;
		this.variantByOrdinal[ordinal] = variantIndex;
	}
	preparePartitionJobs() {
		if (this.partitionJobsCache) return this.partitionJobsCache;
		for (let ordinal = 0; ordinal < this.records.length; ordinal += 1) if (!this.processed[ordinal]) {
			if (this.options.unresolvedPolicy === "retain") this.standalone.push(ordinal);
			this.unresolvedRecords += 1;
		}
		const partitions = /* @__PURE__ */ new Map();
		this.variants.forEach((variant, index) => {
			const values = partitions.get(variant.partition);
			if (values) values.push(index);
			else partitions.set(variant.partition, [index]);
		});
		this.partitionJobsCache = [...partitions.values()].map((group, id) => ({
			id,
			variantIndices: Int32Array.from(group),
			variants: group.map((index) => ({
				sequence: this.arena.decode(this.variants[index].location),
				representative: this.variants[index].representative,
				count: this.variants[index].count
			})),
			options: this.options
		}));
		return this.partitionJobsCache;
	}
	finish(onProgress) {
		const results = this.preparePartitionJobs().map(runDenoisePartitionJob);
		return this.finishWithPartitionResults(results, onProgress);
	}
	finishWithPartitionResults(results, onProgress, reportVariantProgress = true) {
		const jobs = this.preparePartitionJobs();
		const resultById = new Map(results.map((result) => [result.id, result]));
		let candidateComparisons = 0;
		let truncated = 0;
		let indelMergedVariants = 0;
		let substitutionMergedVariants = 0;
		const variantWork = Math.max(1, jobs.reduce((total, job) => total + (this.options.mode === "fad" ? job.variants.length + job.variants.reduce((count, variant) => count + (variant.count >= this.options.minimumParentCount ? 1 : 0), 0) : job.variants.length), 0));
		let processedVariantWork = 0;
		if (reportVariantProgress) onProgress?.(0, variantWork, "variants");
		for (const job of jobs) {
			const result = resultById.get(job.id);
			if (!result || result.targets.length !== job.variantIndices.length) throw new Error(`Denoising partition ${job.id} did not return a complete result.`);
			for (let localIndex = 0; localIndex < job.variantIndices.length; localIndex += 1) {
				const localTarget = result.targets[localIndex];
				if (localTarget < 0 || localTarget >= job.variantIndices.length) throw new Error(`Denoising partition ${job.id} returned an invalid target.`);
				this.variants[job.variantIndices[localIndex]].target = job.variantIndices[localTarget];
			}
			candidateComparisons += result.comparisons;
			truncated += result.truncated;
			indelMergedVariants += result.indelMergedVariants;
			substitutionMergedVariants += result.substitutionMergedVariants;
			processedVariantWork += result.variantWork;
			if (reportVariantProgress) onProgress?.(Math.min(processedVariantWork, variantWork), variantWork, "variants");
		}
		if (reportVariantProgress && !jobs.length) onProgress?.(variantWork, variantWork, "variants");
		const representatives = new Int32Array(this.records.length);
		representatives.fill(-1);
		const counts = new Uint32Array(this.records.length);
		const finalizeTotal = Math.max(1, this.records.length * 3 + this.variants.length + this.standalone.length);
		let finalized = 0;
		const finalizeStride = Math.max(1, Math.floor(finalizeTotal / 500));
		const reportFinalize = (increment = 1) => {
			finalized += increment;
			if (finalized >= finalizeTotal || finalized % finalizeStride === 0) onProgress?.(Math.min(finalized, finalizeTotal), finalizeTotal, "finalize");
		};
		onProgress?.(0, finalizeTotal, "finalize");
		for (let ordinal = 0; ordinal < this.records.length; ordinal += 1) {
			const variantIndex = this.variantByOrdinal[ordinal];
			if (variantIndex >= 0) {
				const target = this.variants[this.variants[variantIndex].target];
				representatives[ordinal] = target.representative;
			}
			reportFinalize();
		}
		for (const variant of this.variants) {
			const target = this.variants[variant.target];
			counts[target.representative] += variant.count;
			reportFinalize();
		}
		for (const group of this.retainedAmbiguous.values()) {
			counts[group.representative] = group.count;
			group.members.forEach((ordinal) => {
				representatives[ordinal] = group.representative;
			});
		}
		for (const ordinal of this.standalone) {
			counts[ordinal] = Math.max(1, Math.floor(this.records[ordinal].inputCount ?? 1));
			representatives[ordinal] = ordinal;
			reportFinalize();
		}
		let uniqueRecords = 0;
		for (const count of counts) {
			if (count > 0) uniqueRecords += 1;
			reportFinalize();
		}
		const largestGroups = largestCountGroups(counts);
		const warnings = [];
		if (this.excludedAmbiguous) warnings.push(`${this.excludedAmbiguous.toLocaleString()} records containing ambiguous nucleotide symbols were excluded to match the selected policy.`);
		if (this.unresolvedRecords) warnings.push(`${this.unresolvedRecords.toLocaleString()} records without a usable trimmed sequence or both V/J calls were ${this.options.unresolvedPolicy === "retain" ? "retained unchanged" : "discarded from the downstream representative set"}.`);
		if (truncated) warnings.push(`${truncated.toLocaleString()} variants reached the candidate cap; increase it before treating this denoising result as complete.`);
		let inputAbundance = 0;
		for (const record of this.records) {
			inputAbundance += Math.max(1, Math.floor(record.inputCount ?? 1));
			reportFinalize();
		}
		return {
			mode: this.options.mode,
			key: "trimmed",
			algorithm: this.options.mode === "fad" ? `FAD-compatible corrected 6-mer / method ${this.options.fadMethod}` : this.options.mode === "indel" ? "Indel-aware bounded edit model" : "Conservative exact-neighbor error model",
			inputRecords: this.records.length,
			inputAbundance,
			uniqueRecords,
			collapsedRecords: this.records.length - uniqueRecords,
			representatives,
			counts,
			largestGroups,
			partitions: jobs.length,
			candidateComparisons,
			indelMergedVariants,
			substitutionMergedVariants,
			excludedAmbiguous: this.excludedAmbiguous,
			unresolvedRecords: this.unresolvedRecords,
			warnings
		};
	}
	processFadPartition(group, onVariant) {
		const maximumSquared = Math.max(0, Math.floor(12 * this.options.fadNeighborThreshold + 1e-9));
		const blockCount = Math.max(1, maximumSquared + 1);
		const profiles = /* @__PURE__ */ new Map();
		const sequences = /* @__PURE__ */ new Map();
		for (const index of group) {
			const sequence = this.arena.decode(this.variants[index].location);
			sequences.set(index, sequence);
			profiles.set(index, kmerProfile(sequence, blockCount));
		}
		const ordered = [...group].sort((left, right) => this.variants[right].count - this.variants[left].count || this.variants[left].representative - this.variants[right].representative);
		const accepted = [];
		const index = /* @__PURE__ */ new Map();
		const candidates = /* @__PURE__ */ new Set();
		let comparisons = 0;
		let truncated = 0;
		const addAccepted = (variantIndex) => {
			accepted.push(variantIndex);
			profiles.get(variantIndex).hashes.forEach((hash, block) => {
				const key = `${block}:${hash}`;
				const values = index.get(key);
				if (values) values.push(variantIndex);
				else index.set(key, [variantIndex]);
			});
		};
		for (const variantIndex of ordered.filter((value) => this.variants[value].count >= this.options.minimumParentCount)) {
			const profile = profiles.get(variantIndex);
			candidates.clear();
			profile.hashes.forEach((hash, block) => {
				if (candidates.size >= this.options.maxCandidatesPerVariant) return;
				for (const candidate of index.get(`${block}:${hash}`) ?? []) {
					if (candidates.size >= this.options.maxCandidatesPerVariant) break;
					candidates.add(candidate);
				}
			});
			if (candidates.size >= this.options.maxCandidatesPerVariant) truncated += 1;
			const neighbors = [];
			for (const candidate of candidates) {
				comparisons += 1;
				const distance = kmerSquaredDistance(profile, profiles.get(candidate), maximumSquared);
				if (distance <= maximumSquared) neighbors.push({
					index: candidate,
					distance
				});
			}
			if (!neighbors.length) {
				addAccepted(variantIndex);
				onVariant?.();
				continue;
			}
			neighbors.sort((left, right) => this.variants[right.index].count - this.variants[left.index].count || left.distance - right.distance || this.variants[left.index].representative - this.variants[right.index].representative);
			const parent = neighbors[0].index;
			const child = this.variants[variantIndex];
			const lambda = this.variants[parent].count / Math.max(Number.MIN_VALUE, this.options.expectedZeroErrorFraction) * this.options.errorRate;
			const adjusted = Math.min(1, poissonStrictUpperTail(child.count, lambda) * (sequences.get(variantIndex)?.length ?? 1));
			if (this.options.fadMethod === 2 && adjusted < this.options.alpha) addAccepted(variantIndex);
			else child.target = parent;
			onVariant?.();
		}
		if (!accepted.length && ordered.length) addAccepted(ordered[0]);
		const acceptedSet = new Set(accepted);
		const vpTree = buildKmerVpTree(accepted, profiles, () => {
			comparisons += 1;
		});
		for (const variantIndex of ordered) {
			if (acceptedSet.has(variantIndex)) {
				this.variants[variantIndex].target = variantIndex;
				onVariant?.();
				continue;
			}
			this.variants[variantIndex].target = vpTree ? nearestKmerPoint(vpTree, variantIndex, profiles, (point) => this.variants[point].count, () => {
				comparisons += 1;
			}) : variantIndex;
			onVariant?.();
		}
		return {
			comparisons,
			truncated,
			indelMergedVariants: 0,
			substitutionMergedVariants: 0
		};
	}
	processIndelPartition(group, onVariant) {
		const distanceLimit = this.options.maximumEditDistance;
		const blockCount = distanceLimit + 1;
		const sequences = new Map(group.map((index) => [index, this.arena.decode(this.variants[index].location)]));
		const ordered = [...group].sort((left, right) => this.variants[right].count - this.variants[left].count || this.variants[left].representative - this.variants[right].representative);
		const parentIndex = /* @__PURE__ */ new Map();
		const shortParentsByLength = /* @__PURE__ */ new Map();
		const profileEdit = createBoundedEditProfiler(distanceLimit);
		let comparisons = 0;
		let truncated = 0;
		let indelMergedVariants = 0;
		let substitutionMergedVariants = 0;
		const candidates = /* @__PURE__ */ new Set();
		const addParent = (variantIndex) => {
			if (this.variants[variantIndex].count < this.options.minimumParentCount) return;
			const sequence = sequences.get(variantIndex);
			const segments = indexedEditSegments(sequence, blockCount);
			if (segments.some((segment) => segment.length === 0)) {
				const values = shortParentsByLength.get(sequence.length);
				if (values) values.push(variantIndex);
				else shortParentsByLength.set(sequence.length, [variantIndex]);
				return;
			}
			for (const segment of segments) {
				const key = `${sequence.length}:${segment.index}:${segment.value}`;
				const values = parentIndex.get(key);
				if (values) values.push(variantIndex);
				else parentIndex.set(key, [variantIndex]);
			}
		};
		for (const variantIndex of ordered) {
			const child = this.variants[variantIndex];
			const sequence = sequences.get(variantIndex);
			candidates.clear();
			let capped = false;
			const addCandidates = (values) => {
				for (const candidate of values) {
					if (candidates.has(candidate)) continue;
					if (candidates.size >= this.options.maxCandidatesPerVariant) {
						capped = true;
						return;
					}
					candidates.add(candidate);
				}
			};
			const minimumParentLength = Math.max(1, sequence.length - distanceLimit);
			const maximumParentLength = sequence.length + distanceLimit;
			for (let parentLength = minimumParentLength; parentLength <= maximumParentLength && !capped; parentLength += 1) {
				addCandidates(shortParentsByLength.get(parentLength) ?? []);
				for (let segmentIndex = 0; segmentIndex < blockCount && !capped; segmentIndex += 1) {
					const parentStart = Math.floor(segmentIndex * parentLength / blockCount);
					const segmentLength = Math.floor((segmentIndex + 1) * parentLength / blockCount) - parentStart;
					if (!segmentLength) continue;
					const firstStart = Math.max(0, parentStart - distanceLimit);
					const lastStart = Math.min(sequence.length - segmentLength, parentStart + distanceLimit);
					for (let queryStart = firstStart; queryStart <= lastStart && !capped; queryStart += 1) addCandidates(parentIndex.get(`${parentLength}:${segmentIndex}:${sequence.slice(queryStart, queryStart + segmentLength)}`) ?? []);
				}
			}
			if (capped) truncated += 1;
			let best = -1;
			let bestProfile = null;
			for (const candidate of candidates) {
				const parent = this.variants[candidate];
				if (parent.count <= child.count) continue;
				comparisons += 1;
				const profile = profileEdit(sequences.get(candidate), sequence);
				if (!profile || profile.distance < 1) continue;
				const indels = profile.insertions + profile.deletions;
				let plausible = false;
				if (indels > 0) plausible = parent.count / child.count >= this.options.minimumIndelParentRatio;
				else {
					const exactErrorProbability = (this.options.errorRate / 3) ** profile.substitutions * (1 - this.options.errorRate) ** Math.max(0, sequence.length - profile.substitutions);
					const lambda = parent.count * exactErrorProbability;
					plausible = Math.min(1, poissonStrictUpperTail(child.count, lambda) * alternativeCount(sequence.length, profile.substitutions)) >= this.options.alpha;
				}
				if (!plausible) continue;
				const bestParent = best >= 0 ? this.variants[best] : null;
				if (!bestProfile || profile.distance < bestProfile.distance || profile.distance === bestProfile.distance && profile.substitutions < bestProfile.substitutions || profile.distance === bestProfile.distance && profile.substitutions === bestProfile.substitutions && parent.count > (bestParent?.count ?? -1) || profile.distance === bestProfile.distance && profile.substitutions === bestProfile.substitutions && parent.count === bestParent?.count && parent.representative < (bestParent?.representative ?? Number.POSITIVE_INFINITY)) {
					best = candidate;
					bestProfile = profile;
				}
			}
			if (best >= 0 && bestProfile) {
				child.target = best;
				if (bestProfile.insertions + bestProfile.deletions > 0) indelMergedVariants += 1;
				else substitutionMergedVariants += 1;
			} else {
				child.target = variantIndex;
				addParent(variantIndex);
			}
			onVariant?.();
		}
		return {
			comparisons,
			truncated,
			indelMergedVariants,
			substitutionMergedVariants
		};
	}
	processConservativePartition(group, onVariant) {
		const distanceLimit = this.options.maximumHammingDistance;
		const blockCount = distanceLimit + 1;
		const sequences = new Map(group.map((index) => [index, this.arena.decode(this.variants[index].location)]));
		const ordered = [...group].sort((left, right) => this.variants[right].count - this.variants[left].count || this.variants[left].representative - this.variants[right].representative);
		const parentIndex = /* @__PURE__ */ new Map();
		let comparisons = 0;
		let truncated = 0;
		let substitutionMergedVariants = 0;
		const candidates = /* @__PURE__ */ new Set();
		const addParent = (variantIndex) => {
			if (this.variants[variantIndex].count < this.options.minimumParentCount) return;
			const sequence = sequences.get(variantIndex);
			sequenceBlocks(sequence, blockCount).forEach((block, blockIndex) => {
				const key = `${sequence.length}:${blockIndex}:${block}`;
				const values = parentIndex.get(key);
				if (values) values.push(variantIndex);
				else parentIndex.set(key, [variantIndex]);
			});
		};
		for (const variantIndex of ordered) {
			const child = this.variants[variantIndex];
			const sequence = sequences.get(variantIndex);
			candidates.clear();
			sequenceBlocks(sequence, blockCount).forEach((block, blockIndex) => {
				if (candidates.size >= this.options.maxCandidatesPerVariant) return;
				for (const candidate of parentIndex.get(`${sequence.length}:${blockIndex}:${block}`) ?? []) {
					if (candidates.size >= this.options.maxCandidatesPerVariant) break;
					candidates.add(candidate);
				}
			});
			if (candidates.size >= this.options.maxCandidatesPerVariant) truncated += 1;
			let best = -1;
			let bestLambda = -1;
			let bestDistance = Number.POSITIVE_INFINITY;
			for (const candidate of candidates) {
				comparisons += 1;
				const distance = hammingDistanceWithin(sequence, sequences.get(candidate), distanceLimit, false);
				if (distance < 1 || distance > distanceLimit) continue;
				const lambda = this.variants[candidate].count * ((this.options.errorRate / 3) ** distance * (1 - this.options.errorRate) ** Math.max(0, sequence.length - distance));
				if (Math.min(1, poissonStrictUpperTail(child.count, lambda) * alternativeCount(sequence.length, distance)) >= this.options.alpha && (lambda > bestLambda || lambda === bestLambda && distance < bestDistance)) {
					best = candidate;
					bestLambda = lambda;
					bestDistance = distance;
				}
			}
			if (best >= 0) {
				child.target = best;
				substitutionMergedVariants += 1;
			} else {
				child.target = variantIndex;
				addParent(variantIndex);
			}
			onVariant?.();
		}
		return {
			comparisons,
			truncated,
			indelMergedVariants: 0,
			substitutionMergedVariants
		};
	}
};
var UnionFind = class {
	parent;
	rank;
	constructor(size) {
		this.parent = new Int32Array(size);
		this.rank = new Uint8Array(size);
		for (let index = 0; index < size; index += 1) this.parent[index] = index;
	}
	find(value) {
		let root = value;
		while (this.parent[root] !== root) root = this.parent[root];
		while (this.parent[value] !== value) {
			const next = this.parent[value];
			this.parent[value] = root;
			value = next;
		}
		return root;
	}
	union(left, right) {
		let a = this.find(left);
		let b = this.find(right);
		if (a === b) return;
		if (this.rank[a] < this.rank[b]) [a, b] = [b, a];
		this.parent[b] = a;
		if (this.rank[a] === this.rank[b]) this.rank[a] += 1;
	}
};
function blocks(value, count) {
	const result = [];
	for (let index = 0; index < count; index += 1) {
		const start = Math.floor(index * value.length / count);
		const end = Math.floor((index + 1) * value.length / count);
		result.push({
			index,
			value: value.slice(start, end)
		});
	}
	return result;
}
function compatibleRecords(left, right, options) {
	if (datasetScopeKey(left, options.scope ?? "global") !== datasetScopeKey(right, options.scope ?? "global")) return false;
	if (options.requireSameLocus && left.locus !== right.locus) return false;
	return callsCompatible(left.vCall, right.vCall, options.callResolution, options.ambiguity) && callsCompatible(left.jCall, right.jCall, options.callResolution, options.ambiguity);
}
function recordIndexTokens(record, options) {
	const v = callSet(record.vCall, options.callResolution, options.ambiguity);
	const j = callSet(record.jCall, options.callResolution, options.ambiguity);
	if (!v.length || !j.length) return [];
	const scope = `${datasetScopeKey(record, options.scope ?? "global")}\u0001`;
	if (options.ambiguity === "strict") return [`${scope}${v.join("+")}\u0001${j.join("+")}`];
	return v.flatMap((vCall) => j.map((jCall) => `${scope}${vCall}\u0001${jCall}`));
}
function assignLineages(records, options, dedup, activeMask, doubleDMask) {
	const union = new UnionFind(records.length);
	const bucket = /* @__PURE__ */ new Map();
	const exact = /* @__PURE__ */ new Map();
	let candidateComparisons = 0;
	let truncatedCandidates = 0;
	let assignedRecords = 0;
	let activeAbundance = 0;
	for (let index = 0; index < records.length; index += 1) {
		const record = records[index];
		const weight = activeMask && !activeMask[index] ? 0 : dedup ? dedup.counts[index] : Math.max(1, Math.floor(records[index].inputCount ?? 1));
		if (!weight) continue;
		activeAbundance += weight;
		const cdr3 = normalizeNt(record.cdr3Nt);
		const tokens = recordIndexTokens(record, options);
		if (!cdr3 || !tokens.length || options.productiveOnly && !record.productive) continue;
		assignedRecords += weight;
		const distanceLimit = Math.floor((1 - options.identity) * cdr3.length + 1e-9);
		const blockCount = Math.max(1, Math.min(cdr3.length, distanceLimit + 1));
		const prefix = `${options.requireSameLocus ? record.locus : "*"}\u0000${cdr3.length}\u0000`;
		const exactSignature = `${prefix}${tokens.join("|")}\u0000${cdr3}`;
		const exactRepresentative = exact.get(exactSignature);
		if (exactRepresentative !== void 0) {
			union.union(index, exactRepresentative);
			continue;
		}
		exact.set(exactSignature, index);
		const candidates = /* @__PURE__ */ new Set();
		for (const token of tokens) {
			for (const block of blocks(cdr3, blockCount)) {
				const key = `${prefix}${token}\u0000${block.index}\u0000${block.value}`;
				for (const candidate of bucket.get(key) ?? []) {
					candidates.add(candidate);
					if (candidates.size >= options.maxCandidateComparisons) break;
				}
				if (candidates.size >= options.maxCandidateComparisons) break;
			}
			if (candidates.size >= options.maxCandidateComparisons) break;
		}
		if (candidates.size >= options.maxCandidateComparisons) truncatedCandidates += 1;
		for (const candidate of candidates) {
			candidateComparisons += 1;
			const other = records[candidate];
			if (!compatibleRecords(record, other, options)) continue;
			if (hammingDistanceWithin(cdr3, normalizeNt(other.cdr3Nt), distanceLimit, false) <= distanceLimit) union.union(index, candidate);
		}
		for (const token of tokens) for (const block of blocks(cdr3, blockCount)) {
			const key = `${prefix}${token}\u0000${block.index}\u0000${block.value}`;
			const values = bucket.get(key);
			if (values) values.push(index);
			else bucket.set(key, [index]);
		}
	}
	const eligible = new Uint8Array(records.length);
	const uniqueByRoot = new Uint32Array(records.length);
	const abundanceByRoot = new Float64Array(records.length);
	const representativeByRoot = new Int32Array(records.length);
	representativeByRoot.fill(-1);
	for (let index = 0; index < records.length; index += 1) {
		const weight = activeMask && !activeMask[index] ? 0 : dedup ? dedup.counts[index] : Math.max(1, Math.floor(records[index].inputCount ?? 1));
		if (!weight) continue;
		const record = records[index];
		if (!record.cdr3Nt || !recordIndexTokens(record, options).length || options.productiveOnly && !record.productive) continue;
		const root = union.find(index);
		eligible[index] = 1;
		uniqueByRoot[root] += 1;
		abundanceByRoot[root] += weight;
		if (representativeByRoot[root] < 0) representativeByRoot[root] = index;
	}
	const assignments = new Int32Array(records.length);
	const rootIds = new Int32Array(records.length);
	const top = [];
	const summaryLimit = 1e4;
	const sizeBins = /* @__PURE__ */ new Uint32Array(7);
	const vUsageMap = /* @__PURE__ */ new Map();
	const jUsageMap = /* @__PURE__ */ new Map();
	let lineageCount = 0;
	const smaller = (left, right) => left.abundance < right.abundance || left.abundance === right.abundance && left.uniqueMembers < right.uniqueMembers || left.abundance === right.abundance && left.uniqueMembers === right.uniqueMembers && left.id > right.id;
	const retainTop = (summary) => {
		if (top.length < summaryLimit) {
			top.push(summary);
			let index = top.length - 1;
			while (index > 0) {
				const parent = Math.floor((index - 1) / 2);
				if (!smaller(top[index], top[parent])) break;
				[top[index], top[parent]] = [top[parent], top[index]];
				index = parent;
			}
			return;
		}
		if (!smaller(top[0], summary)) return;
		top[0] = summary;
		let index = 0;
		while (true) {
			const left = index * 2 + 1;
			const right = left + 1;
			let smallest = index;
			if (left < top.length && smaller(top[left], top[smallest])) smallest = left;
			if (right < top.length && smaller(top[right], top[smallest])) smallest = right;
			if (smallest === index) break;
			[top[index], top[smallest]] = [top[smallest], top[index]];
			index = smallest;
		}
	};
	const addUsage = (map, call, abundance) => {
		const value = map.get(call);
		if (value) {
			value.lineages += 1;
			value.abundance += abundance;
		} else map.set(call, {
			lineages: 1,
			abundance
		});
	};
	for (let root = 0; root < records.length; root += 1) {
		if (!uniqueByRoot[root]) continue;
		const id = ++lineageCount;
		rootIds[root] = id;
		const representativeOrdinal = representativeByRoot[root];
		const representative = records[representativeOrdinal];
		const abundance = abundanceByRoot[root];
		const bin = abundance === 1 ? 0 : abundance <= 3 ? 1 : abundance <= 9 ? 2 : abundance <= 24 ? 3 : abundance <= 99 ? 4 : abundance <= 499 ? 5 : 6;
		sizeBins[bin] += 1;
		const vCalls = callSet(representative.vCall, options.callResolution, options.ambiguity);
		const jCalls = callSet(representative.jCall, options.callResolution, options.ambiguity);
		vCalls.forEach((call) => addUsage(vUsageMap, call, abundance));
		jCalls.forEach((call) => addUsage(jUsageMap, call, abundance));
		retainTop({
			id,
			representativeOrdinal,
			uniqueMembers: uniqueByRoot[root],
			abundance,
			locus: representative.locus,
			vCalls,
			jCalls,
			cdr3Length: representative.cdr3Nt.length,
			representativeCdr3Nt: representative.cdr3Nt,
			representativeCdr3Aa: representative.cdr3Aa,
			studyScope: options.scope ?? "global",
			studyGroup: datasetScopeValue(representative, options.scope ?? "global"),
			sampleIds: [],
			sampleCounts: [],
			subjectIds: [],
			timepoints: [],
			compartments: [],
			doubleDPositiveMembers: 0,
			doubleDPositiveAbundance: 0
		});
	}
	for (let index = 0; index < records.length; index += 1) if (eligible[index]) assignments[index] = rootIds[union.find(index)];
	if (dedup) for (let index = 0; index < dedup.representatives.length; index += 1) {
		const representative = dedup.representatives[index];
		assignments[index] = representative >= 0 ? assignments[representative] : 0;
	}
	const summaryById = new Map(top.map((summary) => [summary.id, summary]));
	const samplesById = /* @__PURE__ */ new Map();
	const sampleCountsById = /* @__PURE__ */ new Map();
	const subjectsById = /* @__PURE__ */ new Map();
	const timepointsById = /* @__PURE__ */ new Map();
	const compartmentsById = /* @__PURE__ */ new Map();
	const addValue = (map, id, value) => {
		if (!value) return;
		const values = map.get(id);
		if (values) values.add(value);
		else map.set(id, /* @__PURE__ */ new Set([value]));
	};
	for (let index = 0; index < records.length; index += 1) {
		const representative = dedup ? dedup.representatives[index] : index;
		if (representative < 0 || activeMask && !activeMask[representative]) continue;
		const sourceWeight = Math.max(1, Math.floor(records[index].inputCount ?? 1));
		const id = assignments[representative];
		if (!id || !summaryById.has(id)) continue;
		const record = records[index];
		addValue(samplesById, id, record.sampleId);
		if (record.sampleId) {
			let bySample = sampleCountsById.get(id);
			if (!bySample) {
				bySample = /* @__PURE__ */ new Map();
				sampleCountsById.set(id, bySample);
			}
			const counts = bySample.get(record.sampleId) ?? {
				representatives: /* @__PURE__ */ new Set(),
				abundance: 0
			};
			counts.representatives.add(representative);
			counts.abundance += sourceWeight;
			bySample.set(record.sampleId, counts);
		}
		addValue(subjectsById, id, record.subjectId);
		addValue(timepointsById, id, record.timepoint);
		addValue(compartmentsById, id, record.compartment);
	}
	for (let index = 0; index < records.length; index += 1) {
		const activeWeight = dedup ? dedup.counts[index] : Math.max(1, Math.floor(records[index].inputCount ?? 1));
		if (!activeWeight || activeMask && !activeMask[index] || !doubleDMask?.[index]) continue;
		const id = assignments[index];
		if (summaryById.has(id)) {
			const summary = summaryById.get(id);
			summary.doubleDPositiveMembers = (summary.doubleDPositiveMembers ?? 0) + 1;
			summary.doubleDPositiveAbundance = (summary.doubleDPositiveAbundance ?? 0) + activeWeight;
		}
	}
	for (const summary of top) {
		summary.sampleIds = [...samplesById.get(summary.id) ?? []].sort((a, b) => a.localeCompare(b, void 0, { numeric: true }));
		summary.sampleCounts = [...sampleCountsById.get(summary.id) ?? []].map(([sampleId, counts]) => ({
			sampleId,
			uniqueMembers: counts.representatives.size,
			abundance: counts.abundance
		})).sort((left, right) => right.abundance - left.abundance || left.sampleId.localeCompare(right.sampleId, void 0, { numeric: true }));
		summary.subjectIds = [...subjectsById.get(summary.id) ?? []].sort((a, b) => a.localeCompare(b, void 0, { numeric: true }));
		summary.timepoints = [...timepointsById.get(summary.id) ?? []].sort((a, b) => a.localeCompare(b, void 0, { numeric: true }));
		summary.compartments = [...compartmentsById.get(summary.id) ?? []].sort((a, b) => a.localeCompare(b, void 0, { numeric: true }));
	}
	const usage = (map) => [...map.entries()].map(([call, value]) => ({
		call,
		...value
	})).sort((a, b) => b.abundance - a.abundance || b.lineages - a.lineages || a.call.localeCompare(b.call));
	return {
		assignments,
		summaries: top.sort((a, b) => b.abundance - a.abundance || b.uniqueMembers - a.uniqueMembers || a.id - b.id),
		lineageCount,
		sizeHistogram: [
			"1",
			"2–3",
			"4–9",
			"10–24",
			"25–99",
			"100–499",
			"500+"
		].map((label, index) => ({
			label,
			count: sizeBins[index]
		})),
		vUsage: usage(vUsageMap),
		jUsage: usage(jUsageMap),
		assignedRecords,
		unassignedRecords: activeAbundance - assignedRecords,
		candidateComparisons,
		truncatedCandidates
	};
}
function parseFasta(text, aligned = false) {
	const records = [];
	let name = "";
	let sequence = [];
	const commit = () => {
		if (!name) return;
		const joined = sequence.join("").toUpperCase().replaceAll(".", "-").replace(/\s/g, "");
		const normalized = aligned ? joined.replace(/[^ACGTNRYKMSWBDHV-]/g, "N") : normalizeNt(joined);
		records.push({
			name,
			sequence: normalized
		});
	};
	for (const line of text.split(/\r?\n/)) if (line.startsWith(">")) {
		commit();
		name = line.slice(1).trim().split(/\s+/, 1)[0] ?? "";
		sequence = [];
	} else if (name) sequence.push(line);
	commit();
	return records;
}
function nameAliases(name) {
	const gene = name.split("|").filter(Boolean).find((field) => /(?:IGH|IGK|IGL|TRA|TRB|TRD|TRG)[VJ]/i.test(field));
	return [...new Set([
		name,
		gene ?? "",
		name.split(/\s+/, 1)[0] ?? ""
	].filter(Boolean))];
}
function prepareReferenceMsa(text) {
	const records = parseFasta(text, true);
	if (records.length < 2) throw new Error("CHMMAIRRa requires at least two aligned reference sequences.");
	const length = records[0].sequence.length;
	if (!length || records.some((record) => record.sequence.length !== length)) throw new Error("Reference MSA records must all have the same aligned length.");
	const lookup = /* @__PURE__ */ new Map();
	records.forEach((record, index) => nameAliases(record.name).forEach((alias) => lookup.set(alias, index)));
	return {
		names: records.map((record) => record.name),
		sequences: records.map((record) => record.sequence),
		lookup,
		length
	};
}
function topCallName(value) {
	return value.split(",", 1)[0]?.trim() ?? "";
}
function threadSequenceToMsa(sequenceAlignment, germlineAlignment, call, msa) {
	const callName = topCallName(call);
	let referenceIndex = msa.lookup.get(callName);
	if (referenceIndex === void 0) {
		const gene = callName.replace(/\*.*$/, "");
		referenceIndex = [...msa.lookup.entries()].find(([name]) => name.replace(/\*.*$/, "") === gene)?.[1];
	}
	if (referenceIndex === void 0) throw new Error(`Reference MSA is missing ${callName || "the assigned allele"}.`);
	const gappedFullReference = msa.sequences[referenceIndex];
	const degappedReference = gappedFullReference.replaceAll("-", "");
	const localGermline = germlineAlignment.toUpperCase().replaceAll(".", "-");
	const localQuery = sequenceAlignment.toUpperCase().replaceAll(".", "-").replace(/[^ACGT-]/g, "N");
	const degappedLocal = localGermline.replaceAll("-", "");
	const matchStart = degappedReference.indexOf(degappedLocal);
	if (matchStart < 0) throw new Error(`The ${callName} local germline alignment cannot be located in its MSA reference.`);
	const withoutInsertions = [];
	for (let index = 0; index < localGermline.length; index += 1) if (localGermline[index] !== "-") withoutInsertions.push(localQuery[index] ?? "N");
	const query = [..."N".repeat(matchStart), ...withoutInsertions];
	const threaded = new Array(msa.length).fill("-");
	let position = 0;
	for (let index = 0; index < msa.length; index += 1) {
		if (position >= query.length) break;
		if (gappedFullReference[index] !== "-") threaded[index] = query[position++] ?? "N";
	}
	return threaded.join("");
}
function chmmairraDistanceFromReference(query, germline) {
	const normalizedQuery = query.toUpperCase().replaceAll(".", "-");
	const normalizedGermline = germline.toUpperCase().replaceAll(".", "-");
	let mismatches = 0;
	for (let index = 0; index < Math.min(normalizedQuery.length, normalizedGermline.length); index += 1) if (normalizedQuery[index] !== normalizedGermline[index]) mismatches += 1;
	return mismatches;
}
function emission(reference, observation, mutation) {
	if (observation === "-" || observation === "N") return 1;
	return reference === observation ? 1 - mutation : mutation / 3;
}
function approximateEmissions(msa, observation, mutations) {
	const values = new Float64Array(msa.sequences.length * msa.length);
	for (let state = 0; state < msa.sequences.length; state += 1) for (let site = 0; site < msa.length; site += 1) values[state * msa.length + site] = emission(msa.sequences[state][site], observation[site], mutations[state]);
	return values;
}
function estimateMutationProbabilities(msa, observation, prior, baseMutation) {
	const states = msa.sequences.length;
	const length = msa.length;
	const mutations = new Float64Array(states);
	mutations.fill(baseMutation);
	if (states < 2 || length < 1) return mutations;
	const emissions = approximateEmissions(msa, observation, mutations);
	const alpha = new Float64Array(states * length);
	const beta = new Float64Array(states * length);
	const scale = new Float64Array(length);
	scale[0] = 1;
	const switchProbability = prior / length;
	const same = 1 - switchProbability;
	const different = switchProbability / (states - 1);
	for (let state = 0; state < states; state += 1) alpha[state * length] = emissions[state * length] / states;
	for (let site = 0; site < length - 1; site += 1) {
		let total = 0;
		for (let state = 0; state < states; state += 1) total += alpha[state * length + site];
		let nextTotal = 0;
		for (let state = 0; state < states; state += 1) {
			const value = ((total - alpha[state * length + site]) * different + alpha[state * length + site] * same) * emissions[state * length + site + 1];
			alpha[state * length + site + 1] = value;
			nextTotal += value;
		}
		scale[site + 1] = nextTotal ? 1 / nextTotal : 1;
		for (let state = 0; state < states; state += 1) alpha[state * length + site + 1] *= scale[site + 1];
	}
	for (let state = 0; state < states; state += 1) beta[state * length + length - 1] = scale[length - 1];
	for (let site = length - 2; site >= 0; site -= 1) {
		let total = 0;
		for (let state = 0; state < states; state += 1) total += beta[state * length + site + 1] * emissions[state * length + site + 1];
		for (let state = 0; state < states; state += 1) {
			const own = beta[state * length + site + 1] * emissions[state * length + site + 1];
			beta[state * length + site] = (different * (total - own) + same * own) * scale[site];
		}
	}
	const mutated = new Float64Array(states);
	const sameCount = new Float64Array(states);
	mutated.fill(2);
	sameCount.fill(10);
	for (let site = 0; site < length; site += 1) {
		if (observation[site] === "-" || observation[site] === "N") continue;
		let normalization = 0;
		for (let state = 0; state < states; state += 1) normalization += alpha[state * length + site] * beta[state * length + site];
		if (!normalization) continue;
		for (let state = 0; state < states; state += 1) {
			const posterior = alpha[state * length + site] * beta[state * length + site] / normalization;
			if (observation[site] === msa.sequences[state][site]) sameCount[state] += posterior;
			else mutated[state] += posterior;
		}
	}
	for (let state = 0; state < states; state += 1) mutations[state] = mutated[state] / (mutated[state] + sameCount[state]);
	return mutations;
}
function approximateProbability(msa, observation, prior, mutations) {
	const states = msa.sequences.length;
	if (states < 2) return 0;
	const emissions = approximateEmissions(msa, observation, mutations);
	const nonChimeric = new Float64Array(states);
	const chimeric = new Float64Array(states);
	const switchProbability = prior / msa.length;
	const same = 1 - switchProbability;
	const different = switchProbability / (states - 1);
	for (let state = 0; state < states; state += 1) nonChimeric[state] = emissions[state * msa.length] / states;
	for (let site = 0; site < msa.length - 1; site += 1) {
		let total = 0;
		for (let state = 0; state < states; state += 1) total += nonChimeric[state] + chimeric[state];
		let nextTotal = 0;
		for (let state = 0; state < states; state += 1) {
			const emissionValue = emissions[state * msa.length + site + 1];
			const own = nonChimeric[state] + chimeric[state];
			chimeric[state] = ((total - own) * different + chimeric[state] * same) * emissionValue;
			nonChimeric[state] = nonChimeric[state] * same * emissionValue;
			nextTotal += chimeric[state] + nonChimeric[state];
		}
		const scaling = nextTotal ? 1 / nextTotal : 1;
		for (let state = 0; state < states; state += 1) {
			chimeric[state] *= scaling;
			nonChimeric[state] *= scaling;
		}
	}
	let chimericTotal = 0;
	let total = 0;
	for (let state = 0; state < states; state += 1) {
		chimericTotal += chimeric[state];
		total += chimeric[state] + nonChimeric[state];
	}
	return total ? chimericTotal / total : 0;
}
function fullProbability(msa, observation, options) {
	const references = msa.sequences.length;
	const rates = options.mutationRates.length ? options.mutationRates : [.005];
	const rateCount = rates.length;
	const states = references * rateCount;
	if (references < 2) return 0;
	const nonChimeric = new Float64Array(states);
	const chimeric = new Float64Array(states);
	const nextNonChimeric = new Float64Array(states);
	const nextChimeric = new Float64Array(states);
	const switchProbability = options.priorProbability / msa.length;
	const self = 1 - switchProbability - options.mutationSwitchProbability;
	const differentReference = switchProbability / ((references - 1) * rateCount);
	const differentMutation = rateCount === 1 ? 0 : options.mutationSwitchProbability / (rateCount - 1);
	const emissionAt = (state, site) => {
		const reference = Math.floor(state / rateCount);
		const rate = state % rateCount;
		return emission(msa.sequences[reference][site], observation[site], rates[rate]);
	};
	for (let state = 0; state < states; state += 1) nonChimeric[state] = emissionAt(state, 0) / states;
	for (let site = 0; site < msa.length - 1; site += 1) {
		let total = 0;
		for (let state = 0; state < states; state += 1) total += nonChimeric[state] + chimeric[state];
		let nextTotal = 0;
		for (let reference = 0; reference < references; reference += 1) {
			let nonReference = 0;
			let chimReference = 0;
			for (let rate = 0; rate < rateCount; rate += 1) {
				const state = reference * rateCount + rate;
				nonReference += nonChimeric[state];
				chimReference += chimeric[state];
			}
			for (let rate = 0; rate < rateCount; rate += 1) {
				const state = reference * rateCount + rate;
				const emit = emissionAt(state, site + 1);
				nextNonChimeric[state] = (nonChimeric[state] * self + (nonReference - nonChimeric[state]) * differentMutation) * emit;
				nextChimeric[state] = (chimeric[state] * self + (chimReference - chimeric[state]) * differentMutation + (total - chimReference - nonReference) * differentReference) * emit;
				nextTotal += nextNonChimeric[state] + nextChimeric[state];
			}
		}
		const scaling = nextTotal ? 1 / nextTotal : 1;
		for (let state = 0; state < states; state += 1) {
			nonChimeric[state] = nextNonChimeric[state] * scaling;
			chimeric[state] = nextChimeric[state] * scaling;
		}
	}
	let chimericTotal = 0;
	let total = 0;
	for (let state = 0; state < states; state += 1) {
		chimericTotal += chimeric[state];
		total += chimeric[state] + nonChimeric[state];
	}
	return total ? chimericTotal / total : 0;
}
function fullViterbi(msa, observation, options) {
	const references = msa.sequences.length;
	const rates = options.mutationRates.length ? options.mutationRates : [.005];
	const rateCount = rates.length;
	const states = references * rateCount;
	const length = msa.length;
	if (references < 2) return {
		startingReference: msa.names[0] ?? "",
		recombinations: [],
		referencePath: options.tracePath ? new Int32Array(length) : void 0
	};
	const previous = new Float64Array(states);
	const next = new Float64Array(states);
	const from = new Int32Array(states * length);
	const switchProbability = options.priorProbability / length;
	const logSelf = Math.log(1 - switchProbability - options.mutationSwitchProbability);
	const logDifferentReference = Math.log(switchProbability / ((references - 1) * rateCount));
	const logDifferentMutation = rateCount === 1 ? Number.NEGATIVE_INFINITY : Math.log(options.mutationSwitchProbability / (rateCount - 1));
	const emissionAt = (state, site) => emission(msa.sequences[Math.floor(state / rateCount)][site], observation[site], rates[state % rateCount]);
	for (let state = 0; state < states; state += 1) {
		previous[state] = Math.log(emissionAt(state, 0) / states);
		from[state * length] = state;
	}
	for (let site = 0; site < length - 1; site += 1) {
		const bestInReference = new Int32Array(references);
		for (let reference = 0; reference < references; reference += 1) {
			let best = reference * rateCount;
			for (let rate = 1; rate < rateCount; rate += 1) {
				const state = reference * rateCount + rate;
				if (previous[state] > previous[best]) best = state;
			}
			bestInReference[reference] = best;
		}
		let bestReference = 0;
		let secondReference = references > 1 ? 1 : 0;
		if (previous[bestInReference[secondReference]] > previous[bestInReference[bestReference]]) [bestReference, secondReference] = [secondReference, bestReference];
		for (let reference = 2; reference < references; reference += 1) if (previous[bestInReference[reference]] > previous[bestInReference[bestReference]]) {
			secondReference = bestReference;
			bestReference = reference;
		} else if (previous[bestInReference[reference]] > previous[bestInReference[secondReference]]) secondReference = reference;
		for (let reference = 0; reference < references; reference += 1) {
			const outside = bestInReference[reference === bestReference ? secondReference : bestReference];
			for (let rate = 0; rate < rateCount; rate += 1) {
				const state = reference * rateCount + rate;
				let selected = state;
				let score = previous[state] + logSelf;
				const referenceScore = previous[outside] + logDifferentReference;
				if (referenceScore > score) {
					selected = outside;
					score = referenceScore;
				}
				if (rateCount > 1) {
					let otherRate = reference * rateCount + (rate === 0 ? 1 : 0);
					for (let candidateRate = 0; candidateRate < rateCount; candidateRate += 1) {
						const candidate = reference * rateCount + candidateRate;
						if (candidate !== state && previous[candidate] > previous[otherRate]) otherRate = candidate;
					}
					const mutationScore = previous[otherRate] + logDifferentMutation;
					if (mutationScore > score) {
						selected = otherRate;
						score = mutationScore;
					}
				}
				from[state * length + site + 1] = selected;
				next[state] = score + Math.log(emissionAt(state, site + 1));
			}
		}
		previous.set(next);
	}
	let current = 0;
	for (let state = 1; state < states; state += 1) if (previous[state] > previous[current]) current = state;
	const recombinations = [];
	const referencePath = options.tracePath ? new Int32Array(length) : void 0;
	if (referencePath) referencePath[length - 1] = Math.floor(current / rateCount);
	for (let site = length - 1; site >= 1; site -= 1) {
		const parent = from[current * length + site];
		const leftReference = Math.floor(parent / rateCount);
		const rightReference = Math.floor(current / rateCount);
		if (leftReference !== rightReference) recombinations.push({
			position: site + 1,
			left: msa.names[leftReference],
			right: msa.names[rightReference]
		});
		current = parent;
		if (referencePath) referencePath[site - 1] = Math.floor(current / rateCount);
	}
	return {
		startingReference: msa.names[Math.floor(current / rateCount)],
		recombinations: recombinations.reverse(),
		referencePath
	};
}
function approximateViterbi(msa, observation, prior, mutations, tracePath = false) {
	const states = msa.sequences.length;
	const length = msa.length;
	const emissions = approximateEmissions(msa, observation, mutations);
	const previous = new Float64Array(states);
	const next = new Float64Array(states);
	const from = new Int32Array(states * length);
	const same = Math.log(1 - prior / length);
	const different = Math.log(prior / length / (states - 1));
	for (let state = 0; state < states; state += 1) {
		previous[state] = Math.log(emissions[state * length] / states);
		from[state * length] = state;
	}
	for (let site = 0; site < length - 1; site += 1) {
		let best = 0;
		for (let state = 1; state < states; state += 1) if (previous[state] > previous[best]) best = state;
		for (let state = 0; state < states; state += 1) {
			const stay = previous[state] + same;
			const change = previous[best] + different;
			const remain = state === best || stay > change;
			from[state * length + site + 1] = remain ? state : best;
			next[state] = (remain ? stay : change) + Math.log(emissions[state * length + site + 1]);
		}
		previous.set(next);
	}
	let current = 0;
	for (let state = 1; state < states; state += 1) if (previous[state] > previous[current]) current = state;
	const reverse = [];
	const referencePath = tracePath ? new Int32Array(length) : void 0;
	if (referencePath) referencePath[length - 1] = current;
	for (let site = length - 1; site >= 1; site -= 1) {
		const parent = from[current * length + site];
		if (parent !== current) reverse.push({
			position: site + 1,
			left: msa.names[parent],
			right: msa.names[current]
		});
		current = parent;
		if (referencePath) referencePath[site - 1] = current;
	}
	return {
		startingReference: msa.names[current],
		recombinations: reverse.reverse(),
		referencePath
	};
}
function runChmm(msa, threadedObservation, localSequenceAlignment, localGermlineAlignment, options) {
	if (threadedObservation.length !== msa.length) throw new Error("Threaded query and reference MSA lengths differ.");
	const dfr = chmmairraDistanceFromReference(localSequenceAlignment, localGermlineAlignment);
	if (options.method === "DB") {
		const path = options.detailed ? fullViterbi(msa, threadedObservation, options) : {
			startingReference: "",
			recombinations: []
		};
		return {
			probability: fullProbability(msa, threadedObservation, options),
			dfr,
			...path
		};
	}
	const mutations = estimateMutationProbabilities(msa, threadedObservation, options.priorProbability, options.baseMutationProbability);
	return {
		probability: approximateProbability(msa, threadedObservation, options.priorProbability, mutations),
		dfr,
		...options.detailed ? approximateViterbi(msa, threadedObservation, options.priorProbability, mutations, options.tracePath) : {
			startingReference: "",
			recombinations: []
		}
	};
}
//#endregion
export { annotateDoubleDBatch as _, deduplicate as a, parseFasta as c, runDenoisePartitionJob as d, runExactDedupJob as f, annotateAirrBatch as g, DEFAULT_PIPELINE_PLAN as h, createExactDedupPlan as i, prepareReferenceMsa as l, threadSequenceToMsa as m, assignLineages as n, finishExactDedupPlan as o, sequenceFingerprint as p, chmmairraDistanceFromReference as r, normalizeNt as s, DenoiseAccumulator as t, runChmm as u, datasetScopeValue as v, stableDatasetSeed as y };
