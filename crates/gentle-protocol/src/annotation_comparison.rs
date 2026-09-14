//! Source-neutral transcript-start annotation comparison contracts.

use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;

pub const TRANSCRIPT_START_ANNOTATION_COMPARISON_SCHEMA: &str =
    "gentle.transcript_start_annotation_comparison.v1";
pub const TRANSCRIPT_START_ANNOTATION_COMPARISON_RECEIPT_SCHEMA: &str =
    "gentle.transcript_start_annotation_comparison_receipt.v1";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TranscriptAnnotationSourceIdentity {
    pub genome_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(flatten)]
    pub metadata: BTreeMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptStartNearestCluster {
    pub tss_1based: usize,
    pub delta_bp_transcript_direction: isize,
    pub absolute_delta_bp: usize,
    pub relationship: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct PrimaryTranscriptStartCluster {
    pub promoter_id: String,
    pub tss_1based: usize,
    pub transcript_ids: Vec<String>,
    pub transcript_count: usize,
    pub selected_for_reporter: bool,
    pub nearest_secondary: TranscriptStartNearestCluster,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SecondaryTranscriptStartCluster {
    pub tss_1based: usize,
    pub transcript_ids: Vec<String>,
    pub transcript_count: usize,
    pub nearest_primary: TranscriptStartNearestCluster,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct GeneTranscriptStartAnnotationComparison {
    pub gene_symbol: String,
    pub chromosome: String,
    pub strand: String,
    pub primary_transcript_count: usize,
    pub secondary_transcript_count: usize,
    pub primary_tss_cluster_count: usize,
    pub secondary_tss_cluster_count: usize,
    pub exact_shared_tss_count: usize,
    pub primary_clusters: Vec<PrimaryTranscriptStartCluster>,
    pub secondary_clusters: Vec<SecondaryTranscriptStartCluster>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct AnnotationComparisonFileBinding {
    pub path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sha1: Option<String>,
    pub sha256: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptStartAnnotationComparisonInputs {
    pub primary_tss_report: AnnotationComparisonFileBinding,
    pub catalog: AnnotationComparisonFileBinding,
    pub secondary_manifest: AnnotationComparisonFileBinding,
    pub secondary_annotation: AnnotationComparisonFileBinding,
    pub secondary_transcript_index: AnnotationComparisonFileBinding,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TranscriptStartAnnotationComparisonReport {
    pub schema: String,
    pub coordinate_system: String,
    pub comparison_rule: String,
    pub primary: TranscriptAnnotationSourceIdentity,
    pub secondary: TranscriptAnnotationSourceIdentity,
    pub inputs: TranscriptStartAnnotationComparisonInputs,
    pub genes: Vec<GeneTranscriptStartAnnotationComparison>,
    pub producer: TranscriptStartAnnotationComparisonProducer,
    pub non_claims: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptStartAnnotationComparisonPageBinding {
    pub gene_symbol: String,
    pub path: String,
    pub sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptStartAnnotationComparisonProducer {
    pub revision: String,
    pub script: String,
    pub script_sha256: String,
    pub gentle_cli: String,
    pub gentle_cli_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptStartAnnotationComparisonReceipt {
    pub schema: String,
    pub report_schema: String,
    pub report_sha256: String,
    pub pages: Vec<TranscriptStartAnnotationComparisonPageBinding>,
    pub outputs: BTreeMap<String, String>,
    pub producer: TranscriptStartAnnotationComparisonProducer,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_delta_round_trips_without_becoming_an_unsigned_distance() {
        let value = serde_json::json!({
            "tss_1based": 2000,
            "delta_bp_transcript_direction": -10,
            "absolute_delta_bp": 10,
            "relationship": "nearest"
        });
        let parsed: TranscriptStartNearestCluster = serde_json::from_value(value).unwrap();
        assert_eq!(parsed.delta_bp_transcript_direction, -10);
        assert_eq!(parsed.absolute_delta_bp, 10);
    }
}
