use std_shims::{vec, String};
use thiserror::Error;

use crate::core::air::{Component, Components};
use crate::core::channel::{Channel, MerkleChannel};
use crate::core::circle::CirclePoint;
use crate::core::fields::qm31::{SecureField, SECURE_EXTENSION_DEGREE};
use crate::core::fri::FriVerificationError;
use crate::core::pcs::CommitmentSchemeVerifier;
use crate::core::proof::StarkProof;
use crate::core::vcs::verifier::MerkleVerificationError;
pub const PREPROCESSED_TRACE_IDX: usize = 0;

pub fn verify<MC: MerkleChannel>(
    components: &[&dyn Component],
    channel: &mut MC::C,
    commitment_scheme: &mut CommitmentSchemeVerifier<MC>,
    proof: StarkProof<MC::H>,
) -> Result<(), VerificationError> {
    let total_start = std::time::Instant::now();
    
    // Setup phase
    let phase_timer = std::time::Instant::now();
    let n_preprocessed_columns = commitment_scheme.trees[PREPROCESSED_TRACE_IDX]
        .column_log_sizes
        .len();

    let components = Components {
        components: components.to_vec(),
        n_preprocessed_columns,
    };
    
    let composition_log_degree = components.composition_log_degree_bound();
    log::info!(
        "Setup: Composition polynomial log degree bound: {}, duration: {:?}",
        composition_log_degree,
        phase_timer.elapsed()
    );
    
    let random_coeff = channel.draw_secure_felt();

    // Commitment verification phase
    let phase_timer = std::time::Instant::now();
    commitment_scheme.commit(
        *proof.commitments.last().unwrap(),
        &[composition_log_degree; SECURE_EXTENSION_DEGREE],
        channel,
    );
    log::info!("Commitment verification took: {:?}", phase_timer.elapsed());

    // OODS point generation
    let phase_timer = std::time::Instant::now();
    let oods_point = CirclePoint::<SecureField>::get_random_point(channel);
    log::info!("OODS point generation took: {:?}", phase_timer.elapsed());

    // Sample points computation
    let phase_timer = std::time::Instant::now();
    let mut sample_points = components.mask_points(oods_point);
    sample_points.push(vec![vec![oods_point]; SECURE_EXTENSION_DEGREE]);

    let sample_points_by_column = sample_points.as_cols_ref().flatten();
    let n_columns = sample_points_by_column.len();
    let total_sample_points = sample_points_by_column.into_iter().flatten().count();
    
    log::info!("Sampling {} columns with {} total sample points, took: {:?}", 
        n_columns, total_sample_points, phase_timer.elapsed());

    // OODS evaluation and verification
    let phase_timer = std::time::Instant::now();
    let composition_oods_eval =
        proof
            .extract_composition_oods_eval()
            .ok_or(VerificationError::InvalidStructure(
                std_shims::ToString::to_string(&"Unexpected sampled_values structure"),
            ))?;

    let expected_composition_eval = components.eval_composition_polynomial_at_point(
        oods_point,
        &proof.sampled_values,
        random_coeff,
    );

    if composition_oods_eval != expected_composition_eval {
        log::error!("OODS verification failed - values don't match");
        return Err(VerificationError::OodsNotMatching);
    }
    
    log::info!("OODS verification took: {:?}", phase_timer.elapsed());

    // Final verification
    let phase_timer = std::time::Instant::now();
    let result = commitment_scheme.verify_values(sample_points, proof.0, channel);
    log::info!("Final verification took: {:?}, success: {}", phase_timer.elapsed(), result.is_ok());
    
    log::info!("Total verification time: {:?}", total_start.elapsed());
    
    result
}

#[derive(Clone, Debug, Error)]
pub enum VerificationError {
    #[error("Proof has invalid structure: {0}.")]
    InvalidStructure(String),
    #[error(transparent)]
    Merkle(#[from] MerkleVerificationError),
    #[error(
        "The composition polynomial OODS value does not match the trace OODS values
    (DEEP-ALI failure)."
    )]
    OodsNotMatching,
    #[error(transparent)]
    Fri(#[from] FriVerificationError),
    #[error("Proof of work verification failed.")]
    ProofOfWork,
}
