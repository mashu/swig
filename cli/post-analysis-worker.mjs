import { m as runExactDedupJob, p as runDenoisePartitionJob } from "./chunks/post-analysis-core-LY5d8K8i.mjs";
import { parentPort } from "node:worker_threads";
//#region cli-src/post-analysis-worker.mjs
const send = (message, transfer) => parentPort ? parentPort.postMessage(message, transfer) : globalThis.postMessage(message, transfer);
const receive = (message) => {
	try {
		if (message.kind === "denoise") {
			const result = runDenoisePartitionJob(message.job);
			send({
				id: message.id,
				result
			}, [result.targets.buffer]);
		} else if (message.kind === "exact") {
			const result = runExactDedupJob(message.job);
			send({
				id: message.id,
				result
			}, [
				result.representatives.buffer,
				result.representativeOrdinals.buffer,
				result.representativeCounts.buffer
			]);
		} else throw new Error(`Unknown post-analysis worker task: ${String(message.kind)}`);
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
