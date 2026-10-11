import { readFileSync } from "node:fs";
import { parentPort } from "node:worker_threads";
import { AlignerEngine } from "../src/aligner-engine.ts";

const encoder=new TextEncoder();
const decoder=new TextDecoder();
let engine=null;
let doubleD=false;

async function initialize(message){
  // The pool compiles the module once and shares it; compile here only if it did not.
  const module=message.module??await WebAssembly.compile(readFileSync(message.wasmPath));
  const created=await AlignerEngine.create(module,message.settings,message.references);
  engine=created.engine;
  doubleD=Boolean(message.settings.doubleD&&message.settings.doubleD.mode!=="off");
  return {records:created.records,summary:created.engine.summary};
}

function annotate(message){
  if(!engine)throw new Error("The aligner worker is not initialized.");
  const result=engine.annotate(encoder.encode(message.text),message.count,doubleD);
  return {
    header:result.header,
    body:decoder.decode(result.body),
    count:result.count,
    doubleDHeader:result.doubleDHeader??"",
    doubleDBody:result.doubleDBody?decoder.decode(result.doubleDBody):"",
  };
}

const send=(message)=>parentPort?parentPort.postMessage(message):globalThis.postMessage(message);
const receive=async(message)=>{
  try{
    const result=message.type==="init"?await initialize(message):annotate(message);
    send({id:message.id,...result});
  }catch(error){send({id:message.id,error:error instanceof Error?error.message:String(error)});}
};

if(parentPort)parentPort.on("message",receive);
else globalThis.addEventListener("message",(event)=>void receive(event.data));
