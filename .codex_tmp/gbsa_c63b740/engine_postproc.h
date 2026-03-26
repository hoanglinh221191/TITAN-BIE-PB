#ifndef ENGINE_POSTPROC_H
#define ENGINE_POSTPROC_H

// engine_postproc.h — Post-loop analysis orchestration
// Extracted from impl_engine.cpp for compilation speed and clarity.
// Each function owns its own timing and logs to the shared fstream.

#include "arguments.h"
#include "cas_analysis.h"
#include "dccm_network.h"
#include "info.h"
#include "ml_extractor.h"
#include "nucleic_profiler.h"
#include "titan_kinetics.h"
#include "titan_llps.h"
#include "titan_msm.h"
#include "titan_pca.h"
#include "water_analysis.h"

#include <chrono>
#include <fstream>
#include <memory>

#include "gbsa_intrin.h" // GBSA_Aggregator complete type

// ---- 1. GBSA + ML + CAS + Nucleic Acid post-processing -------------------
void run_gbsa_postproc(const arguments &args, size_t tot_snap_allfiles,
                       MLFeatureExtractor &ml_extractor,
                       GBSA_Aggregator &global_gbsa_agg,
                       CASAnalyzer &cas_analyzer, NucleicProfiler &na_profiler,
                       bool global_has_nucleic_acid, info &info_output,
                       std::fstream &log);

// ---- 2. TICA + MSM post-processing ---------------------------------------
void run_tica_msm_postproc(const arguments &args,
                           std::unique_ptr<TitanKinetics> &tica_analyzer,
                           std::unique_ptr<TitanMSM> &msm_analyzer,
                           size_t num_threads, info &info_output,
                           std::fstream &log);

// ---- 3. LLPS post-processing ----------------------------------------------
void run_llps_postproc(const arguments &args,
                       std::unique_ptr<TitanLLPS> &llps_analyzer,
                       size_t num_threads, info &info_output,
                       std::fstream &log);

// ---- 4. DCCM + PCA post-processing ----------------------------------------
void run_dccm_pca_postproc(const arguments &args, size_t tot_snap_allfiles,
                           DCCMNetwork &dccm_analyzer, size_t num_threads,
                           std::chrono::milliseconds &time_network,
                           info &info_output, std::fstream &log);

// ---- 5. LSP-MD post-processing ---------------------------------------------
void run_lspmd_postproc(const arguments &args, DCCMNetwork &dccm_analyzer,
                        size_t num_threads, info &info_output,
                        std::fstream &log);

// ---- 6. DHAM post-processing -----------------------------------------------
#include "titan_dham.h"
void run_dham_postproc(const arguments &args,
                       std::unique_ptr<TitanMSM> &msm_analyzer,
                       size_t num_threads, info &info_output,
                       std::fstream &log);

#endif // ENGINE_POSTPROC_H
