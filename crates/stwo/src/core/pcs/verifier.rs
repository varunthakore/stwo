use core::iter::zip;

use itertools::Itertools;
use std_shims::Vec;
use log::info;  // Add this import

use super::super::circle::CirclePoint;
use super::super::fields::qm31::SecureField;
use super::super::fri::{CirclePolyDegreeBound, FriVerifier};
use super::quotients::{fri_answers, PointSample};
use super::utils::TreeVec;
use super::PcsConfig;
use crate::core::channel::{Channel, MerkleChannel};
use crate::core::pcs::quotients::CommitmentSchemeProof;
use crate::core::vcs::verifier::MerkleVerifier;
use crate::core::vcs::MerkleHasher;
use crate::core::verifier::VerificationError;
use crate::core::ColumnVec;

/// The verifier side of a FRI polynomial commitment scheme. See [super].
#[derive(Default)]
pub struct CommitmentSchemeVerifier<MC: MerkleChannel> {
    pub trees: TreeVec<MerkleVerifier<MC::H>>,
    pub config: PcsConfig,
}

impl<MC: MerkleChannel> CommitmentSchemeVerifier<MC> {
    pub fn new(config: PcsConfig) -> Self {
        Self {
            trees: TreeVec::default(),
            config,
        }
    }

    /// A [TreeVec<ColumnVec>] of the log sizes of each column in each commitment tree.
    fn column_log_sizes(&self) -> TreeVec<ColumnVec<u32>> {
        self.trees
            .as_ref()
            .map(|tree| tree.column_log_sizes.clone())
    }

    /// Reads a commitment from the prover.
    pub fn commit(
        &mut self,
        commitment: <MC::H as MerkleHasher>::Hash,
        log_sizes: &[u32],
        channel: &mut MC::C,
    ) {
        MC::mix_root(channel, commitment);
        let extended_log_sizes = log_sizes
            .iter()
            .map(|&log_size| log_size + self.config.fri_config.log_blowup_factor)
            .collect();
        let verifier = MerkleVerifier::new(commitment, extended_log_sizes);
        self.trees.push(verifier);
    }

    pub fn verify_values(
        &self,
        sampled_points: TreeVec<ColumnVec<Vec<CirclePoint<SecureField>>>>,
        proof: CommitmentSchemeProof<MC::H>,
        channel: &mut MC::C,
    ) -> Result<(), VerificationError> {
        let total_start = std::time::Instant::now();
        info!("CommitmentSchemeVerifier::verify_values - Starting verification");
        
        // Channel mixing and random coefficient generation
        let phase_timer = std::time::Instant::now();
        channel.mix_felts(&proof.sampled_values.clone().flatten_cols());
        let random_coeff = channel.draw_secure_felt();
        info!("CommitmentSchemeVerifier::verify_values - Channel mixing and random coeff generation took: {:?}", phase_timer.elapsed());

        // Bounds computation
        let phase_timer = std::time::Instant::now();
        let bounds = self
            .column_log_sizes()
            .flatten()
            .into_iter()
            .sorted()
            .rev()
            .dedup()
            .map(|log_size| {
                CirclePolyDegreeBound::new(log_size - self.config.fri_config.log_blowup_factor)
            })
            .collect_vec();
        info!("CommitmentSchemeVerifier::verify_values - Bounds computation took: {:?}, found {} bounds", phase_timer.elapsed(), bounds.len());

        // FRI commitment phase on OODS quotients
        let phase_timer = std::time::Instant::now();
        let mut fri_verifier =
            FriVerifier::<MC>::commit(channel, self.config.fri_config, proof.fri_proof, bounds)?;
        info!("CommitmentSchemeVerifier::verify_values - FRI commitment phase took: {:?}", phase_timer.elapsed());

        // Verify proof of work
        let phase_timer = std::time::Instant::now();
        channel.mix_u64(proof.proof_of_work);
        if channel.trailing_zeros() < self.config.pow_bits {
            info!("CommitmentSchemeVerifier::verify_values - Proof of work verification failed");
            return Err(VerificationError::ProofOfWork);
        }
        info!("CommitmentSchemeVerifier::verify_values - Proof of work verification took: {:?}", phase_timer.elapsed());

        // Get FRI query positions
        let phase_timer = std::time::Instant::now();
        let query_positions_per_log_size = fri_verifier.sample_query_positions(channel);
        info!("CommitmentSchemeVerifier::verify_values - FRI query positions sampling took: {:?}", phase_timer.elapsed());

        // Verify merkle decommitments
        let phase_timer = std::time::Instant::now();
        self.trees
            .as_ref()
            .zip_eq(proof.decommitments)
            .zip_eq(proof.queried_values.clone())
            .map(|((tree, decommitment), queried_values)| {
                tree.verify(&query_positions_per_log_size, queried_values, decommitment)
            })
            .0
            .into_iter()
            .collect::<Result<(), _>>()?;
        info!("CommitmentSchemeVerifier::verify_values - Merkle decommitments verification took: {:?}", phase_timer.elapsed());

        // Answer FRI queries preparation
        let phase_timer = std::time::Instant::now();
        let samples = sampled_points.zip_cols(proof.sampled_values).map_cols(
            |(sampled_points, sampled_values)| {
                zip(sampled_points, sampled_values)
                    .map(|(point, value)| PointSample { point, value })
                    .collect_vec()
            },
        );

        let n_columns_per_log_size = self.trees.as_ref().map(|tree| &tree.n_columns_per_log_size);

        let fri_answers = fri_answers(
            self.column_log_sizes(),
            samples,
            random_coeff,
            &query_positions_per_log_size,
            proof.queried_values,
            n_columns_per_log_size,
        )?;
        info!("CommitmentSchemeVerifier::verify_values - FRI answers preparation took: {:?}", phase_timer.elapsed());

        // FRI decommitment
        let phase_timer = std::time::Instant::now();
        fri_verifier.decommit(fri_answers)?;
        info!("CommitmentSchemeVerifier::verify_values - FRI decommitment took: {:?}", phase_timer.elapsed());

        info!("CommitmentSchemeVerifier::verify_values - Total verification time: {:?}", total_start.elapsed());
        Ok(())
    }
}
