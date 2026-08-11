#!/usr/bin/env node

import { createHash } from "node:crypto";
import { mkdir, readFile, writeFile } from "node:fs/promises";
import { spawn } from "node:child_process";
import path from "node:path";
import { fileURLToPath } from "node:url";

export const sha256 = (bytes) => `sha256:${createHash("sha256").update(bytes).digest("hex")}`;

export function alignedPrimerEvidence(columns, threePrimeWindowBp) {
  if (columns.length < 16) throw new Error(`expected 16 BLAST columns, observed ${columns.length}`);
  const qstart = Number(columns[6]);
  const qend = Number(columns[7]);
  const qlen = Number(columns[13]);
  const qseq = columns[14].toUpperCase();
  const sseq = columns[15].toUpperCase();
  let queryPos = qstart - 1;
  let mismatches = qstart - 1;
  const observed = new Map();
  for (let i = 0; i < qseq.length; i += 1) {
    if (qseq[i] !== "-") queryPos += 1;
    if (qseq[i] !== "-") {
      const matches = sseq[i] !== "-" && qseq[i] === sseq[i];
      observed.set(queryPos, matches);
      if (!matches) mismatches += 1;
    } else if (sseq[i] !== "-") {
      mismatches += 1;
    }
  }
  mismatches += Math.max(0, qlen - qend);
  let threePrimeMismatches = 0;
  for (let position = Math.max(1, qlen - threePrimeWindowBp + 1); position <= qlen; position += 1) {
    if (observed.get(position) !== true) threePrimeMismatches += 1;
  }
  return {
    query_coverage_fraction: Number(columns[12]) / 100,
    total_mismatches: mismatches,
    three_prime_mismatches: threePrimeMismatches,
    subject_start: Number(columns[8]),
    subject_end: Number(columns[9]),
    orientation: Number(columns[8]) <= Number(columns[9]) ? "plus" : "minus",
  };
}

export function viableHits(tsv, policy) {
  if (!tsv.trim()) return [];
  return tsv.trim().split(/\r?\n/).flatMap((line) => {
    const columns = line.split("\t");
    const evidence = alignedPrimerEvidence(columns, policy.three_prime_window_bp);
    if (evidence.query_coverage_fraction < policy.min_primer_coverage_fraction) return [];
    if (evidence.three_prime_mismatches > policy.max_3prime_mismatches) return [];
    return [{ subject_id: columns[1], ...evidence }];
  });
}

export function collectProducts(forwardHits, reverseHits, maxAmpliconBp) {
  const bySubject = new Map();
  for (const hit of forwardHits) {
    const row = bySubject.get(hit.subject_id) ?? { forward: [], reverse: [] };
    row.forward.push(hit);
    bySubject.set(hit.subject_id, row);
  }
  for (const hit of reverseHits) {
    const row = bySubject.get(hit.subject_id) ?? { forward: [], reverse: [] };
    row.reverse.push(hit);
    bySubject.set(hit.subject_id, row);
  }
  const products = [];
  for (const [subjectId, row] of bySubject) {
    for (const forward of row.forward) for (const reverse of row.reverse) {
      let start;
      let end;
      if (forward.orientation === "plus" && reverse.orientation === "minus") {
        start = Math.min(forward.subject_start, forward.subject_end);
        end = Math.max(reverse.subject_start, reverse.subject_end);
      } else if (forward.orientation === "minus" && reverse.orientation === "plus") {
        start = Math.min(reverse.subject_start, reverse.subject_end);
        end = Math.max(forward.subject_start, forward.subject_end);
      } else continue;
      if (end < start) continue;
      const ampliconLengthBp = end - start + 1;
      if (ampliconLengthBp > maxAmpliconBp) continue;
      products.push({
        subject_id: subjectId,
        start_1based: start,
        end_1based: end,
        amplicon_length_bp: ampliconLengthBp,
        total_primer_mismatches: forward.total_mismatches + reverse.total_mismatches,
        forward,
        reverse,
      });
    }
  }
  return products;
}

const readJson = async (file) => JSON.parse(await readFile(file, "utf8"));

async function run(program, args) {
  return await new Promise((resolve, reject) => {
    const child = spawn(program, args, { stdio: ["ignore", "pipe", "pipe"] });
    const stdout = [];
    const stderr = [];
    child.stdout.on("data", (chunk) => stdout.push(chunk));
    child.stderr.on("data", (chunk) => stderr.push(chunk));
    child.on("error", reject);
    child.on("close", (code, signal) => resolve({
      code,
      signal,
      stdout: Buffer.concat(stdout).toString("utf8"),
      stderr: Buffer.concat(stderr).toString("utf8"),
    }));
  });
}

function intendedSubjects(handoff) {
  return new Set((handoff.intended_target?.expected_products ?? [])
    .filter((row) => row.target_space === "transcriptome_cdna")
    .map((row) => row.subject_id.split(".")[0]));
}

function expandHandoff(document) {
  if (document.schema === "gentle.primer_specificity_handoff.v1") return [document];
  if (document.schema === "gentle.transcript_assay_panel_specificity_handoff.v1") {
    return document.assays.map((row) => row.handoff);
  }
  throw new Error(`unsupported handoff schema ${document.schema}`);
}

async function assessHandoff(handoff, handoffPath, outputDir) {
  const resolvedHandoffPath = handoff.handoff_path ?? handoffPath;
  const handoffBytes = await readFile(resolvedHandoffPath);
  const onDiskHandoff = JSON.parse(handoffBytes.toString("utf8"));
  if (JSON.stringify(onDiskHandoff) !== JSON.stringify(handoff)) {
    throw new Error(`embedded handoff differs from ${resolvedHandoffPath}`);
  }
  const policy = handoff.policy;
  const enrichedOutputs = [];
  for (const command of handoff.commands) {
    if (command.program !== "blastn") throw new Error(`refusing undeclared program ${command.program}`);
    const args = [...command.args];
    const outputIndex = args.indexOf("-out");
    const formatIndex = args.indexOf("-outfmt");
    if (outputIndex < 0 || formatIndex < 0) throw new Error(`${command.command_id} lacks -out or -outfmt`);
    const outputPath = path.join(outputDir, `${handoff.handoff_id}.${command.role}.aligned.tsv`);
    args[outputIndex + 1] = outputPath;
    args[formatIndex + 1] = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen qseq sseq";
    const result = await run(command.program, args);
    if (result.code !== 0) throw new Error(`${command.command_id} failed (${result.code}): ${result.stderr}`);
    const bytes = await readFile(outputPath);
    enrichedOutputs.push({
      role: command.role,
      output_path: outputPath,
      output_size_bytes: bytes.length,
      output_sha256: sha256(bytes),
      text: bytes.toString("utf8"),
    });
  }
  const forwardHits = viableHits(enrichedOutputs.find((row) => row.role === "forward").text, policy);
  const reverseHits = viableHits(enrichedOutputs.find((row) => row.role === "reverse").text, policy);
  const maxAmpliconBp = policy.readiness_max_target_amplicon_bp ?? policy.max_target_amplicon_bp;
  const products = collectProducts(forwardHits, reverseHits, maxAmpliconBp);
  const intended = intendedSubjects(handoff);
  const observedIntended = new Set(products
    .filter((row) => intended.has(row.subject_id.split(".")[0]))
    .map((row) => row.subject_id.split(".")[0]));
  const missingIntended = [...intended].filter((subject) => !observedIntended.has(subject)).sort();
  const unintended = products.filter((row) =>
    !intended.has(row.subject_id.split(".")[0])
    && row.total_primer_mismatches < policy.min_total_mismatches_to_unintended_target);
  const status = missingIntended.length === 0 && unintended.length === 0 ? "pass" : "specificity_fail";
  return {
    schema: "gentle.primer_specificity_linear_fallback_assessment.v1",
    handoff_path: resolvedHandoffPath,
    handoff_sha256: sha256(handoffBytes),
    handoff_id: handoff.handoff_id,
    database_content_fingerprint: handoff.blast_database.content_fingerprint,
    intended_subject_ids: [...intended].sort(),
    observed_intended_subject_ids: [...observedIntended].sort(),
    missing_intended_subject_ids: missingIntended,
    viable_forward_hit_count: forwardHits.length,
    viable_reverse_hit_count: reverseHits.length,
    candidate_product_count: products.length,
    disqualifying_unintended_products: unintended,
    status,
    accepted: status === "pass",
    aligned_blast_outputs: enrichedOutputs.map(({ text: _text, ...row }) => row),
  };
}

async function main() {
  const [requestPath] = process.argv.slice(2);
  if (!requestPath) throw new Error("usage: primer_specificity_linear_fallback.mjs REQUEST.json");
  const requestBytes = await readFile(requestPath);
  const request = JSON.parse(requestBytes.toString("utf8"));
  if (request.schema !== "gentle.primer_specificity_linear_fallback_request.v1") {
    throw new Error(`unsupported request schema ${request.schema}`);
  }
  if (!Array.isArray(request.handoff_paths) || request.handoff_paths.length === 0) {
    throw new Error("request must declare at least one handoff path");
  }
  await mkdir(request.output_dir, { recursive: false });
  const assessments = [];
  for (const handoffPath of request.handoff_paths) {
    const document = await readJson(handoffPath);
    for (const handoff of expandHandoff(document)) {
      assessments.push(await assessHandoff(handoff, handoffPath, request.output_dir));
    }
  }
  const report = {
    schema: "gentle.primer_specificity_linear_fallback_batch.v1",
    request_path: requestPath,
    request_sha256: sha256(requestBytes),
    method: "blast_aligned_strings_linear_subject_pairing_v1",
    generated_utc: new Date().toISOString(),
    assessment_count: assessments.length,
    pass_count: assessments.filter((row) => row.status === "pass").length,
    specificity_fail_count: assessments.filter((row) => row.status === "specificity_fail").length,
    assessments,
  };
  const reportPath = path.join(request.output_dir, "primer_specificity_linear_fallback.batch.json");
  await writeFile(reportPath, `${JSON.stringify(report, null, 2)}\n`, { flag: "wx" });
  console.log(JSON.stringify({ report_path: reportPath, pass_count: report.pass_count, specificity_fail_count: report.specificity_fail_count }));
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main().catch((error) => {
    console.error(error.stack ?? String(error));
    process.exitCode = 1;
  });
}
