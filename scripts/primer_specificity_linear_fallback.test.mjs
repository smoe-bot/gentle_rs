import assert from "node:assert/strict";
import test from "node:test";

import {
  alignedPrimerEvidence,
  collectProducts,
  viableHits,
} from "./primer_specificity_linear_fallback.mjs";

function row({
  subject = "ENST1",
  qstart = 1,
  qend = 20,
  sstart = 100,
  send = 119,
  qcovs = 100,
  qlen = 20,
  qseq = "ACGTACGTACGTACGTACGT",
  sseq = qseq,
} = {}) {
  return [
    "primer", subject, "100", String(qseq.replaceAll("-", "").length), "0", "0",
    String(qstart), String(qend), String(sstart), String(send), "1e-5", "40",
    String(qcovs), String(qlen), qseq, sseq,
  ];
}

const policy = {
  min_primer_coverage_fraction: 0.8,
  max_3prime_mismatches: 0,
  three_prime_window_bp: 5,
};

test("complete aligned strings retain exact terminal evidence", () => {
  const evidence = alignedPrimerEvidence(row(), 5);
  assert.equal(evidence.query_coverage_fraction, 1);
  assert.equal(evidence.total_mismatches, 0);
  assert.equal(evidence.three_prime_mismatches, 0);
  assert.equal(evidence.orientation, "plus");
});

test("unaligned 5-prime bases count as mismatches without invalidating an exact 3-prime window", () => {
  const columns = row({
    qstart: 3,
    qend: 20,
    qcovs: 90,
    qseq: "GTACGTACGTACGTACGT",
    sseq: "GTACGTACGTACGTACGT",
  });
  const [hit] = viableHits(`${columns.join("\t")}\n`, policy);
  assert.equal(hit.total_mismatches, 2);
  assert.equal(hit.three_prime_mismatches, 0);
});

test("a missing or mismatched 3-prime base rejects the hit", () => {
  const columns = row({
    qend: 19,
    qcovs: 95,
    qseq: "ACGTACGTACGTACGTACG",
    sseq: "ACGTACGTACGTACGTACG",
  });
  assert.deepEqual(viableHits(`${columns.join("\t")}\n`, policy), []);
});

test("linear subject grouping emits only inward-facing bounded products", () => {
  const forward = [{
    subject_id: "ENST1", subject_start: 100, subject_end: 119,
    orientation: "plus", total_mismatches: 0,
  }];
  const reverse = [
    {
      subject_id: "ENST1", subject_start: 250, subject_end: 231,
      orientation: "minus", total_mismatches: 1,
    },
    {
      subject_id: "ENST2", subject_start: 250, subject_end: 231,
      orientation: "minus", total_mismatches: 0,
    },
  ];
  const products = collectProducts(forward, reverse, 500);
  assert.equal(products.length, 1);
  assert.equal(products[0].subject_id, "ENST1");
  assert.equal(products[0].amplicon_length_bp, 151);
  assert.equal(products[0].total_primer_mismatches, 1);
});
