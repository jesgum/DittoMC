///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   Ditto.cxx
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Implementation of the Ditto generator.
///

#include "Ditto.h"

#include "DittoSignal.h"
#include "DittoTune.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TParticle.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace Ditto
{

struct TuneCDF {
  std::vector<double> lowEdges;
  std::vector<double> highEdges;
  std::vector<double> cumulative;
};

struct TuneAliasSampler {
  std::vector<double> lowEdges;
  std::vector<double> widths;
  std::vector<double> probability;
  std::vector<std::uint32_t> alias;
};

struct TuneDiscreteAliasSampler {
  std::vector<int> values;
  std::vector<double> probability;
  std::vector<std::uint32_t> alias;
};

struct TuneCompositionTemplateRange {
  std::uint64_t first = 0;
  std::uint64_t last = 0;
};

struct TuneRuntimeSpeciesData {
  int pdg = 0;
  bool charged = false;
  double mass = 0.0;
  double mass2 = 0.0;

  std::vector<TuneAliasSampler> ptGivenActivity;
  std::vector<TuneAliasSampler> etaGivenActivity;
  std::vector<TuneAliasSampler> etaCentralGivenActivity;
  std::vector<TuneAliasSampler> etaOutsideCentralGivenActivity;

  // Exact-Nch mean-pT correction factors learned per species. The correction
  // is split using the same central-charged / other definition as the
  // empirical composition templates.
  std::vector<double> centralPtScaleByNch;
  std::vector<double> otherPtScaleByNch;
};

struct TuneRuntimeData {
  TuneDiscreteAliasSampler nch;

  // Indexed directly by exact Nch.
  std::vector<TuneDiscreteAliasSampler> nSelectedGivenNch;
  std::vector<int> activityClassByNch;

  std::vector<TuneRuntimeSpeciesData> species;

  // Exact (Nch, Nselected) -> empirical composition-template range.
  std::unordered_map<std::uint64_t, TuneCompositionTemplateRange> compositionRanges;

  // Reused scratch buffers to avoid per-event allocations.
  std::vector<int> currentCentralCounts;
  std::vector<int> currentOtherCounts;
};

namespace
{

constexpr double kTwoPi = 6.283185307179586476925286766559;
// Exact-Nch mean-pT bins below this particle count are left uncorrected. The
// raw count histograms are persisted in Tune v6, so this protects the runtime
// model from statistical excursions in very sparse species / multiplicity bins.
constexpr double kMinPtMeanCorrectionParticles = 100.0;

TuneCDF makeCDF(const TH1D& hist,
                double minimum = -std::numeric_limits<double>::infinity(),
                double maximum = std::numeric_limits<double>::infinity())
{
  TuneCDF result;

  double sum = 0.0;
  for (int i = 1; i <= hist.GetNbinsX(); ++i) {
    const double center = hist.GetXaxis()->GetBinCenter(i);
    if (center < minimum || center >= maximum) {
      continue;
    }

    const double weight = std::max(0.0, hist.GetBinContent(i));
    if (weight <= 0.0) {
      continue;
    }

    result.lowEdges.push_back(hist.GetXaxis()->GetBinLowEdge(i));
    result.highEdges.push_back(hist.GetXaxis()->GetBinUpEdge(i));

    sum += weight;
    result.cumulative.push_back(sum);
  }

  if (sum > 0.0) {
    for (auto& value : result.cumulative) {
      value /= sum;
    }
    result.cumulative.back() = 1.0;
  }

  return result;
}

TuneCDF makeYSliceCDF(const TH2D& hist,
                      int xBin,
                      double minimum = -std::numeric_limits<double>::infinity(),
                      double maximum = std::numeric_limits<double>::infinity())
{
  TuneCDF result;

  double sum = 0.0;
  const TAxis* yAxis = hist.GetYaxis();
  for (int i = 1; i <= yAxis->GetNbins(); ++i) {
    const double center = yAxis->GetBinCenter(i);

    if (center < minimum || center >= maximum) {
      continue;
    }

    const double weight = std::max(0.0, hist.GetBinContent(xBin, i));
    if (weight <= 0.0) {
      continue;
    }

    result.lowEdges.push_back(yAxis->GetBinLowEdge(i));
    result.highEdges.push_back(yAxis->GetBinUpEdge(i));

    sum += weight;
    result.cumulative.push_back(sum);
  }

  if (sum > 0.0) {
    for (auto& value : result.cumulative) {
      value /= sum;
    }
    result.cumulative.back() = 1.0;
  }

  return result;
}

TuneCDF makeYSliceCDF(const THnSparseD& hist,
                      int xBin,
                      double minimum = -std::numeric_limits<double>::infinity(),
                      double maximum = std::numeric_limits<double>::infinity())
{
  TuneCDF result;

  double sum = 0.0;
  const TAxis* yAxis = hist.GetAxis(1);
  for (int i = 1; i <= yAxis->GetNbins(); ++i) {
    const double center = yAxis->GetBinCenter(i);

    if (center < minimum || center >= maximum) {
      continue;
    }

    const int coordinate[2] = {xBin, i};
    const int64_t bin = hist.GetBin(coordinate);
    const double weight = std::max(0.0, hist.GetBinContent(bin));
    if (weight <= 0.0) {
      continue;
    }

    result.lowEdges.push_back(yAxis->GetBinLowEdge(i));
    result.highEdges.push_back(yAxis->GetBinUpEdge(i));

    sum += weight;
    result.cumulative.push_back(sum);
  }

  if (sum > 0.0) {
    for (auto& value : result.cumulative) {
      value /= sum;
    }
    result.cumulative.back() = 1.0;
  }

  return result;
}

TuneCDF makeEtaSliceCDF(const TH2D& hist,
                        int xBin,
                        double activityEtaMax,
                        bool central)
{
  TuneCDF result;

  const auto appendInterval = [&](double low,
                                  double high,
                                  double weight,
                                  double& sum) {
    if (high <= low || weight <= 0.0) {
      return;
    }

    result.lowEdges.push_back(low);
    result.highEdges.push_back(high);

    sum += weight;
    result.cumulative.push_back(sum);
  };

  double sum = 0.0;

  const TAxis* yAxis = hist.GetYaxis();
  for (int i = 1; i <= yAxis->GetNbins(); ++i) {
    const double binLow = yAxis->GetBinLowEdge(i);
    const double binHigh = yAxis->GetBinUpEdge(i);
    const double binWidth = binHigh - binLow;
    const double binWeight = std::max(0.0, hist.GetBinContent(xBin, i));

    if (binWeight <= 0.0 || binWidth <= 0.0) {
      continue;
    }

    if (central) {
      const double low = std::max(binLow, -activityEtaMax);
      const double high = std::min(binHigh, activityEtaMax);
      const double fraction = std::max(0.0, high - low) / binWidth;

      appendInterval(low, high, binWeight * fraction, sum);
    } else {
      const double leftLow = binLow;
      const double leftHigh = std::min(binHigh, -activityEtaMax);
      const double leftFraction = std::max(0.0, leftHigh - leftLow) / binWidth;

      appendInterval(leftLow, leftHigh, binWeight * leftFraction, sum);

      const double rightLow = std::max(binLow, activityEtaMax);
      const double rightHigh = binHigh;
      const double rightFraction = std::max(0.0, rightHigh - rightLow) / binWidth;

      appendInterval(rightLow, rightHigh, binWeight * rightFraction, sum);
    }
  }

  if (sum > 0.0) {
    for (auto& value : result.cumulative) {
      value /= sum;
    }
    result.cumulative.back() = 1.0;
  }

  return result;
}

double meanOfCDF(const TuneCDF& cdf)
{
  const std::size_t nBins = cdf.cumulative.size();

  if (nBins == 0 || cdf.lowEdges.size() != nBins || cdf.highEdges.size() != nBins) {
    return 0.0;
  }

  double previous = 0.0;
  double normalization = 0.0;
  double firstMoment = 0.0;

  for (std::size_t i = 0; i < nBins; ++i) {
    const double weight = std::max(0.0, cdf.cumulative[i] - previous);

    previous = cdf.cumulative[i];

    const double center = 0.5 * (cdf.lowEdges[i] + cdf.highEdges[i]);

    normalization += weight;
    firstMoment += weight * center;
  }

  return normalization > 0.0 ? firstMoment / normalization
                             : 0.0;
}

TuneAliasSampler makeAliasSampler(const TuneCDF& cdf)
{
  TuneAliasSampler result;

  const std::size_t nBins = cdf.cumulative.size();

  if (nBins == 0 || cdf.lowEdges.size() != nBins || cdf.highEdges.size() != nBins) {
    return result;
  }

  result.lowEdges = cdf.lowEdges;
  result.widths.resize(nBins);
  result.probability.assign(nBins, 1.0);
  result.alias.resize(nBins);

  std::vector<double> scaled(nBins, 0.0);
  std::vector<std::size_t> small;
  std::vector<std::size_t> large;

  small.reserve(nBins);
  large.reserve(nBins);

  double previous = 0.0;

  for (std::size_t i = 0; i < nBins; ++i) {
    result.widths[i] = cdf.highEdges[i] - cdf.lowEdges[i];

    const double weight = std::max(0.0, cdf.cumulative[i] - previous);

    previous = cdf.cumulative[i];

    scaled[i] = weight * static_cast<double>(nBins);

    result.alias[i] = static_cast<std::uint32_t>(i);

    if (scaled[i] < 1.0) {
      small.push_back(i);
    } else {
      large.push_back(i);
    }
  }

  while (!small.empty() && !large.empty()) {
    const std::size_t smallIndex = small.back();

    small.pop_back();

    const std::size_t largeIndex = large.back();

    large.pop_back();

    result.probability[smallIndex] = std::clamp(scaled[smallIndex], 0.0, 1.0);

    result.alias[smallIndex] = static_cast<std::uint32_t>(largeIndex);

    scaled[largeIndex] = scaled[largeIndex] + scaled[smallIndex] - 1.0;

    if (scaled[largeIndex] < 1.0) {
      small.push_back(largeIndex);
    } else {
      large.push_back(largeIndex);
    }
  }

  for (const std::size_t index : large) {
    result.probability[index] = 1.0;
    result.alias[index] = static_cast<std::uint32_t>(index);
  }

  for (const std::size_t index : small) {
    result.probability[index] = 1.0;
    result.alias[index] = static_cast<std::uint32_t>(index);
  }

  return result;
}

TuneDiscreteAliasSampler makeDiscreteAliasSampler(const TuneCDF& cdf)
{
  TuneDiscreteAliasSampler result;

  const std::size_t nBins = cdf.cumulative.size();

  if (nBins == 0 || cdf.lowEdges.size() != nBins || cdf.highEdges.size() != nBins) {
    return result;
  }

  result.values.resize(nBins);
  result.probability.assign(nBins, 1.0);
  result.alias.resize(nBins);

  std::vector<double> scaled(nBins, 0.0);

  std::vector<std::size_t> small;
  std::vector<std::size_t> large;

  small.reserve(nBins);
  large.reserve(nBins);

  double previous = 0.0;

  for (std::size_t i = 0; i < nBins; ++i) {
    result.values[i] = static_cast<int>(std::llround(0.5 * (cdf.lowEdges[i] + cdf.highEdges[i])));

    const double weight = std::max(0.0, cdf.cumulative[i] - previous);
    previous = cdf.cumulative[i];

    scaled[i] = weight * static_cast<double>(nBins);
    result.alias[i] = static_cast<std::uint32_t>(i);

    if (scaled[i] < 1.0) {
      small.push_back(i);
    } else {
      large.push_back(i);
    }
  }

  while (!small.empty() && !large.empty()) {
    const std::size_t smallIndex = small.back();
    small.pop_back();

    const std::size_t largeIndex = large.back();
    large.pop_back();

    result.probability[smallIndex] = std::clamp(scaled[smallIndex], 0.0, 1.0);
    result.alias[smallIndex] = static_cast<std::uint32_t>(largeIndex);
    scaled[largeIndex] = scaled[largeIndex] + scaled[smallIndex] - 1.0;

    if (scaled[largeIndex] < 1.0) {
      small.push_back(largeIndex);
    } else {
      large.push_back(largeIndex);
    }
  }

  for (const std::size_t index : large) {
    result.probability[index] = 1.0;
    result.alias[index] = static_cast<std::uint32_t>(index);
  }

  for (const std::size_t index : small) {
    result.probability[index] = 1.0;
    result.alias[index] = static_cast<std::uint32_t>(index);
  }

  return result;
}

std::uint64_t compositionKey(int nch, int nSelected)
{
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(nch)) << 32) |
         static_cast<std::uint32_t>(nSelected);
}

} // namespace

Generator::Generator(const Config& config) : mConfig(config),
                                             mRng(config.mSeed)
{
  validateConfig();
  loadTune(mConfig.mTuneFile);

  if (!mConfig.mSignalFile.empty()) {
    constexpr std::uint64_t kSignalSeedOffset = 0x9e3779b97f4a7c15ULL;
    mSignalInjector = new SignalInjector(mConfig.mSignalFile,
                                         mTune->mActivityEdges,
                                         mTune->mFinalStatus,
                                         mConfig.mSeed + kSignalSeedOffset);
  }

  const double meanMultiplicity = mTune->mHNSelected.GetMean();

  const double rmsMultiplicity = mTune->mHNSelected.GetRMS();

  const int initialCapacity = std::max(256, static_cast<int>(std::ceil(meanMultiplicity + 6.0 * rmsMultiplicity)));

  mParticles = new TClonesArray("TParticle", initialCapacity);
}

Generator::~Generator()
{
  if (timingEnabled() && mTiming.mGeneratedEvents > 0) {
    printTimingMetrics(std::cout);
  }

  if (mTree) {
    mOutputFile->cd();
    mTree->Write();
    mOutputFile->Close();
    delete mOutputFile;
    mOutputFile = nullptr;
    mTree = nullptr;
  }

  delete mParticles;
  mParticles = nullptr;

  delete mSignalInjector;
  mSignalInjector = nullptr;

  delete mTuneRuntime;
  mTuneRuntime = nullptr;

  delete mTune;
  mTune = nullptr;
}

const Tune& Generator::generatorCard() const
{
  if (!mTune) {
    throw std::logic_error("Ditto: generator card is not loaded");
  }

  return *mTune;
}

double Generator::averageGenerationTimeUs() const
{
  if (mTiming.mGeneratedEvents == 0) {
    return 0.0;
  }

  return 1.0e6 * mTiming.mGeneration / static_cast<double>(mTiming.mGeneratedEvents);
}

void Generator::resetTimingMetrics()
{
  mTiming = TimingMetrics{};
}

void Generator::printTimingMetrics(std::ostream& os) const
{
  if (mTiming.mGeneratedEvents == 0) {
    os << "Ditto timing: no timed events.\n";
    return;
  }

  const double nEvents = static_cast<double>(mTiming.mGeneratedEvents);
  const double nParticles = static_cast<double>(mTiming.mGeneratedParticles);

  const auto usPerEvent = [nEvents](double seconds) {
    return nEvents > 0.0 ? 1.0e6 * seconds / nEvents : 0.0;
  };
  const auto nsPerParticle = [nParticles](double seconds) {
    return nParticles > 0.0 ? 1.0e9 * seconds / nParticles : 0.0;
  };
  const auto percentOf = [](double part, double total) {
    return total > 0.0 ? 100.0 * part / total : 0.0;
  };

  const double averageUs = usPerEvent(mTiming.mGeneration);
  const double rateHz = averageUs > 0.0 ? 1.0e6 / averageUs : 0.0;
  const double meanMultiplicity = nEvents > 0.0 ? nParticles / nEvents : 0.0;
  const double totalGenerationSec = mTiming.mGeneration;
  const double totalProcessSec = mTiming.mGeneration + mTiming.mTreeFill + mTiming.mPythiaExport;

  os << std::fixed << std::setprecision(3)
     << "\n========== Ditto performance ==========\n"
     << "Events                 : " << mTiming.mGeneratedEvents << "\n"
     << "Particles              : " << mTiming.mGeneratedParticles << "\n"
     << "Mean particles/event   : " << meanMultiplicity << "\n"
     << "Total generation time  : " << totalGenerationSec << " s\n";
  if (mTiming.mTreeFills > 0 || mTiming.mPythiaExports > 0) {
    os << "Total process time     : " << totalProcessSec << " s\n";
  }
  os << "Pure generation        : " << mTiming.mGeneration << " s  (" << percentOf(mTiming.mGeneration, totalProcessSec) << " %)" << "  " << averageUs << " us/event" << "  (" << rateHz << " events/s)\n";

  if (mTune) {
    if (mTune->mNEvents > 0 && mTune->mGenerationTimeSeconds > 0.0) {
      const double tuneAverageUs = 1.0e6 * mTune->mGenerationTimeSeconds / static_cast<double>(mTune->mNEvents);
      const double speedup = averageUs > 0.0 ? tuneAverageUs / averageUs : 0.0;

      os << "Tune teacher generation : "
         << mTune->mGenerationTimeSeconds
         << " s  ("
         << tuneAverageUs
         << " us/event)\n"
         << "Tune speedup            : "
         << speedup
         << "x\n";
    } else {
      os << "Tune teacher generation : timing metadata unavailable\n";
    }
  }

  if (mTiming.mTreeFills > 0) {
    const double treeFillUs = 1.0e6 * mTiming.mTreeFill / static_cast<double>(mTiming.mTreeFills);
    os << "TTree::Fill            : " << mTiming.mTreeFill << " s  (" << percentOf(mTiming.mTreeFill, totalProcessSec) << " %)" << "  " << treeFillUs << " us/fill  (" << percentOf(mTiming.mTreeFill, totalProcessSec) << " %)\n";
  }

  if (mTiming.mPythiaExports > 0) {
    const double pythiaExportUs = 1.0e6 * mTiming.mPythiaExport / static_cast<double>(mTiming.mPythiaExports);
    os << "PYTHIA export          : " << mTiming.mPythiaExport << " s  (" << percentOf(mTiming.mPythiaExport, totalProcessSec) << " %)" << "  " << pythiaExportUs << " us/export  (" << percentOf(mTiming.mPythiaExport, totalProcessSec) << " %)\n";
  }

  if (!detailedTimingEnabled()) {
    os << "Detailed profiling     : disabled\n"
       << "===========================================\n";
    return;
  }

  const double eventAccounted = mTiming.mClearParticles +
                                mTiming.mEventInfo +
                                mTiming.mMultiplicity +
                                mTiming.mExpandArray +
                                mTiming.mParticleLoop;

  const double eventRemainder = std::max(0.0, mTiming.mGeneration - eventAccounted);

  os << "\n-- generation breakdown --\n";

  const auto printEventRow = [&](const char* name, double value) {
    os << std::left << std::setw(24) << name
       << std::right << std::setw(12) << usPerEvent(value) << " us/event"
       << std::setw(10) << percentOf(value, mTiming.mGeneration) << " %\n";
  };

  printEventRow("TClonesArray::Clear", mTiming.mClearParticles);
  printEventRow("event info", mTiming.mEventInfo);
  printEventRow("multiplicity", mTiming.mMultiplicity);
  printEventRow("array expansion", mTiming.mExpandArray);
  printEventRow("particle loop", mTiming.mParticleLoop);
  printEventRow("unaccounted/timer", eventRemainder);

  const double particleAccounted = mTiming.mSpeciesSampling +
                                   mTiming.mConstructedAt +
                                   mTiming.mPtSampling +
                                   mTiming.mEtaSampling +
                                   mTiming.mPhiSampling +
                                   mTiming.mMomentumMath +
                                   mTiming.mParticleSetters;

  const double particleRemainder = std::max(0.0, mTiming.mParticleLoop - particleAccounted);

  os << "\n-- particle-loop breakdown --\n";

  const auto printParticleRow = [&](const char* name, double value) {
    os << std::left << std::setw(24) << name
       << std::right << std::setw(12) << nsPerParticle(value) << " ns/particle"
       << std::setw(10) << percentOf(value, mTiming.mParticleLoop) << " %\n";
  };

  printParticleRow("species sampling", mTiming.mSpeciesSampling);
  printParticleRow("ConstructedAt", mTiming.mConstructedAt);
  printParticleRow("pT sampling", mTiming.mPtSampling);
  printParticleRow("eta sampling", mTiming.mEtaSampling);
  printParticleRow("phi sampling", mTiming.mPhiSampling);
  printParticleRow("momentum math", mTiming.mMomentumMath);
  printParticleRow("TParticle setters", mTiming.mParticleSetters);
  printParticleRow("loop/timer overhead", particleRemainder);

  if (mTiming.mPythiaExports > 0) {
    const double nExports = static_cast<double>(mTiming.mPythiaExports);
    const auto usPerExport = [nExports](double seconds) {
      return nExports > 0.0 ? 1.0e6 * seconds / nExports : 0.0;
    };

    os << "\n-- PYTHIA export breakdown --\n"
       << std::left << std::setw(24) << "event.reset"
       << std::right << std::setw(12) << usPerExport(mTiming.mPythiaReset) << " us/export\n"
       << std::left << std::setw(24) << "append particles"
       << std::right << std::setw(12) << usPerExport(mTiming.mPythiaAppend) << " us/export\n"
       << std::left << std::setw(24) << "system four-vector"
       << std::right << std::setw(12) << usPerExport(mTiming.mPythiaSystemSum) << " us/export\n";
  }

  os << "\nNOTE: detailed profiling inserts several steady_clock::now() calls per\n"
     << "particle. Use enableDetailedTimingMetrics=false for the least intrusive\n"
     << "absolute throughput measurement.\n"
     << "===========================================\n";
}

void Generator::validateConfig() const
{
  if (mConfig.mTuneFile.empty()) {
    throw std::invalid_argument("Ditto: tuneFile is mandatory; the tune is the generator card");
  }
}

void Generator::loadTune(const std::string& fileName)
{
  auto tune = Tune::load(fileName);

  if (tune->mFormatVersion < 8) {
    throw std::runtime_error("Ditto: exact-Nch species-dependent pT correction requires tune formatVersion >= 8");
  }

  if (tune->mAzimuthModel != "uniform") {
    throw std::runtime_error("Ditto: unsupported azimuth model '" + tune->mAzimuthModel + "'");
  }

  auto runtime = std::make_unique<TuneRuntimeData>();

  runtime->nch = makeDiscreteAliasSampler(makeCDF(tune->mPNch,
                                                  tune->mActivityEdges.front(),
                                                  tune->mActivityEdges.back()));

  if (runtime->nch.probability.empty()) {
    throw std::runtime_error("Ditto: tune pNch is empty in the configured activity range");
  }

  // ----------------------------------------------------------------------
  // Exact P(Nselected | Nch), indexed directly by Nch.
  // ----------------------------------------------------------------------

  int maximumNch = 0;

  for (int xBin = 1; xBin <= tune->mPNSelectedGivenNch.GetAxis(0)->GetNbins(); ++xBin) {
    maximumNch = std::max(maximumNch, static_cast<int>(std::llround(tune->mPNSelectedGivenNch.GetAxis(0)->GetBinCenter(xBin))));
  }

  runtime->nSelectedGivenNch.resize(static_cast<std::size_t>(maximumNch + 1));

  runtime->activityClassByNch.assign(static_cast<std::size_t>(maximumNch + 1), -1);

  for (int xBin = 1; xBin <= tune->mPNSelectedGivenNch.GetAxis(0)->GetNbins(); ++xBin) {
    const int nch = static_cast<int>(std::llround(tune->mPNSelectedGivenNch.GetAxis(0)->GetBinCenter(xBin)));

    if (nch < 0) {
      continue;
    }

    const double minimumSelected = std::max<double>(nch, tune->mSelectedMultiplicityEdges.front());

    auto sampler = makeDiscreteAliasSampler(makeYSliceCDF(tune->mPNSelectedGivenNch, xBin, minimumSelected, tune->mSelectedMultiplicityEdges.back()));

    if (sampler.probability.empty() && tune->mPNch.GetBinContent(tune->mPNch.GetXaxis()->FindBin(nch)) > 0.0) {
      throw std::runtime_error("Ditto: empty P(Nselected | Nch=" + std::to_string(nch) +
                               "). Increase selectedMultiplicityEdges and/or "
                               "maxSelectedMultiplicity when tuning.");
    }

    runtime->nSelectedGivenNch[static_cast<std::size_t>(nch)] = std::move(sampler);
    runtime->activityClassByNch[static_cast<std::size_t>(nch)] = tune->activityClass(nch);
  }

  // ----------------------------------------------------------------------
  // Per-species kinematics.
  // ----------------------------------------------------------------------

  const int nActivityClasses = static_cast<int>(tune->mActivityEdges.size()) - 1;

  runtime->species.reserve(tune->numberOfSpecies());

  for (int iSpecies = 0; iSpecies < tune->numberOfSpecies(); ++iSpecies) {
    const auto* entry = tune->speciesAt(iSpecies);

    if (!entry) {
      throw std::runtime_error("Ditto: tune contains an invalid species entry");
    }

    TuneRuntimeSpeciesData speciesData;
    speciesData.pdg = entry->mPdg;
    speciesData.charged = entry->mChargeType != 0;
    speciesData.mass = entry->mMass;
    speciesData.mass2 = entry->mMass * entry->mMass;

    speciesData.ptGivenActivity.reserve(nActivityClasses);
    speciesData.etaGivenActivity.reserve(nActivityClasses);
    speciesData.etaCentralGivenActivity.reserve(nActivityClasses);
    speciesData.etaOutsideCentralGivenActivity.reserve(nActivityClasses);

    std::vector<double> meanPtGivenActivity;
    meanPtGivenActivity.reserve(nActivityClasses);

    for (int i = 1; i <= nActivityClasses; ++i) {
      const auto ptCDF = makeYSliceCDF(entry->mPPtGivenActivity, i);

      meanPtGivenActivity.push_back(meanOfCDF(ptCDF));

      speciesData.ptGivenActivity.push_back(makeAliasSampler(ptCDF));

      speciesData.etaGivenActivity.push_back(makeAliasSampler(makeYSliceCDF(entry->mPEtaGivenActivity, i)));

      speciesData.etaCentralGivenActivity.push_back(makeAliasSampler(makeEtaSliceCDF(entry->mPEtaGivenActivity, i, tune->mActivityEtaMax, true)));

      speciesData.etaOutsideCentralGivenActivity.push_back(makeAliasSampler(makeEtaSliceCDF(entry->mPEtaGivenActivity, i, tune->mActivityEtaMax, false)));
    }

    speciesData.centralPtScaleByNch.assign(static_cast<std::size_t>(maximumNch + 1), 1.0);
    speciesData.otherPtScaleByNch.assign(static_cast<std::size_t>(maximumNch + 1), 1.0);

    const auto fillPtScale = [&](const TH1D& meanPt,
                                 const TH1D& particleCount,
                                 int nch,
                                 double baselineMean,
                                 double& scale) {
      if (baselineMean <= 0.0 || !std::isfinite(baselineMean)) {
        return;
      }

      const int bin = particleCount.GetXaxis()->FindFixBin(nch);

      if (bin < 1 || bin > particleCount.GetNbinsX()) {
        return;
      }

      const double count = particleCount.GetBinContent(bin);

      if (count < kMinPtMeanCorrectionParticles) {
        return;
      }

      const double targetMean = meanPt.GetBinContent(bin);

      if (targetMean <= 0.0 || !std::isfinite(targetMean)) {
        return;
      }

      scale = targetMean / baselineMean;
    };

    for (int nch = 0; nch <= maximumNch; ++nch) {
      const int activityClass = runtime->activityClassByNch[static_cast<std::size_t>(nch)];

      if (activityClass < 0 || activityClass >= nActivityClasses) {
        continue;
      }

      const double baselineMean = meanPtGivenActivity[static_cast<std::size_t>(activityClass)];

      fillPtScale(entry->mHCentralChargedMeanPtVsNch,
                  entry->mHCentralChargedPtCountVsNch,
                  nch,
                  baselineMean,
                  speciesData.centralPtScaleByNch[static_cast<std::size_t>(nch)]);

      fillPtScale(entry->mHOtherMeanPtVsNch,
                  entry->mHOtherPtCountVsNch,
                  nch,
                  baselineMean,
                  speciesData.otherPtScaleByNch[static_cast<std::size_t>(nch)]);
    }

    runtime->species.push_back(std::move(speciesData));
  }

  // ----------------------------------------------------------------------
  // Exact-pair empirical correlated composition templates.
  // ----------------------------------------------------------------------

  runtime->compositionRanges.reserve(static_cast<std::size_t>(tune->numberOfCompositionPairs()));

  for (int iPair = 0; iPair < tune->numberOfCompositionPairs(); ++iPair) {
    const std::uint64_t first = tune->mCompositionPairOffsets[static_cast<std::size_t>(iPair)];

    const std::uint64_t last = tune->mCompositionPairOffsets[static_cast<std::size_t>(iPair) + 1];

    if (last <= first) {
      continue;
    }

    runtime->compositionRanges[compositionKey(tune->mCompositionPairNch[static_cast<std::size_t>(iPair)],
                                              tune->mCompositionPairNSelected[static_cast<std::size_t>(iPair)])] = TuneCompositionTemplateRange{first, last};
  }

  if (runtime->compositionRanges.empty()) {
    throw std::runtime_error("Ditto: tune contains no empirical composition templates");
  }

  runtime->currentCentralCounts.assign(runtime->species.size(), 0);
  runtime->currentOtherCounts.assign(runtime->species.size(), 0);

  if (tune->mCentralChargedCoverageMissingParticles > 0) {
    double totalCentralCharged = 0.0;

    for (int i = 1; i <= tune->mHNch.GetNbinsX(); ++i) {
      totalCentralCharged += tune->mHNch.GetBinContent(i) * tune->mHNch.GetXaxis()->GetBinCenter(i);
    }

    const double missingFraction = totalCentralCharged > 0.0 ? static_cast<double>(tune->mCentralChargedCoverageMissingParticles) / totalCentralCharged : 0.0;

    std::cerr << "Ditto WARNING: the tune species list misses "
              << tune->mCentralChargedCoverageMissingParticles
              << " central charged PYTHIA particles ("
              << 100.0 * missingFraction
              << "% of the central charged sample). "
              << "When such a composition template is selected, Ditto will "
              << "move an existing modeled charged particle from the outer to "
              << "the central component, preserving its species and Nselected.\n";
  }

  delete mTuneRuntime;
  mTuneRuntime = runtime.release();

  delete mTune;
  mTune = tune.release();
}

double Generator::uniform01()
{
  // Convert the upper 53 random bits directly to a double in [0, 1).
  // This avoids constructing std::uniform_real_distribution in hot paths.
  return static_cast<double>(mRng() >> 11) * 0x1.0p-53;
}

double Generator::uniform(double min, double max)
{
  return min + (max - min) * uniform01();
}

EventInfo Generator::makeEventInfo() const
{
  EventInfo info;
  info.mEventNumber = mEventCounter;
  info.mBeamIdA = mTune->mBeamIdA;
  info.mBeamIdB = mTune->mBeamIdB;
  info.mSqrtSNN = mTune->mSqrtSNN;
  return info;
}

int Generator::sampleMultiplicity(EventInfo& info)
{
  if (!mTune || !mTuneRuntime) {
    throw std::logic_error("Ditto: generator card is not loaded");
  }

  const int nch = sampleDiscreteAlias(mTuneRuntime->nch);

  if (nch < 0 || static_cast<std::size_t>(nch) >= mTuneRuntime->nSelectedGivenNch.size()) {
    throw std::runtime_error("Ditto: could not sample a valid Nch from the generator card");
  }

  const int activityClass = mTuneRuntime->activityClassByNch[static_cast<std::size_t>(nch)];

  if (activityClass < 0) {
    throw std::runtime_error("Ditto: sampled Nch is outside the tune activity range");
  }

  const int nSelected = sampleDiscreteAlias(mTuneRuntime->nSelectedGivenNch[static_cast<std::size_t>(nch)]);

  if (nSelected < nch) {
    throw std::runtime_error("Ditto: sampled Nselected is smaller than Nch");
  }

  info.mConditioningNch = nch;
  info.mConditioningNSelected = nSelected;
  info.mActivityClass = activityClass;

  return nSelected;
}

void Generator::sampleComposition(const EventInfo& info,
                                  std::vector<int>& centralCounts,
                                  std::vector<int>& otherCounts)
{
  if (!mTune || !mTuneRuntime) {
    throw std::logic_error("Ditto: no tune runtime cache is loaded");
  }

  const std::uint64_t key = compositionKey(info.mConditioningNch,
                                           info.mConditioningNSelected);

  const auto found = mTuneRuntime->compositionRanges.find(key);

  if (found == mTuneRuntime->compositionRanges.end()) {
    throw std::runtime_error("Ditto: no empirical composition template for sampled (Nch, Nselected) = (" + std::to_string(info.mConditioningNch) + ", " + std::to_string(info.mConditioningNSelected) + ")");
  }

  const auto& range = found->second;

  if (range.last <= range.first) {
    throw std::runtime_error("Ditto: empty empirical composition-template range");
  }

  std::uniform_int_distribution<std::uint64_t> templateDistribution(range.first, range.last - 1);

  const std::uint64_t templateIndex = templateDistribution(mRng);

  const std::size_t nSpecies = mTuneRuntime->species.size();

  const std::uint64_t base = templateIndex * static_cast<std::uint64_t>(nSpecies);

  if (base + nSpecies > mTune->mCompositionCentralCounts.size() || base + nSpecies > mTune->mCompositionOtherCounts.size()) {
    throw std::runtime_error("Ditto: composition-template index is out of range");
  }

  if (centralCounts.size() != nSpecies) {
    centralCounts.resize(nSpecies);
  }

  if (otherCounts.size() != nSpecies) {
    otherCounts.resize(nSpecies);
  }

  int centralTotal = 0;
  int selectedTotal = 0;

  for (std::size_t iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    centralCounts[iSpecies] = static_cast<int>(mTune->mCompositionCentralCounts[static_cast<std::size_t>(base) + iSpecies]);

    otherCounts[iSpecies] = static_cast<int>(mTune->mCompositionOtherCounts[static_cast<std::size_t>(base) + iSpecies]);

    centralTotal += centralCounts[iSpecies];
    selectedTotal += centralCounts[iSpecies] + otherCounts[iSpecies];
  }

  if (selectedTotal != info.mConditioningNSelected) {
    throw std::runtime_error("Ditto: empirical composition template does not reproduce Nselected");
  }

  if (centralTotal > info.mConditioningNch) {
    throw std::runtime_error("Ditto: empirical composition template contains more modeled central charged particles than Nch");
  }

  // The selected species list can miss a tiny fraction of the teacher's
  // central charged particles. Correct only the central/outer assignment:
  // move existing modeled charged particles from the outer component to the
  // central component. This leaves every species total and Nselected exactly
  // unchanged.
  int missingCentral = info.mConditioningNch - centralTotal;

  while (missingCentral > 0) {
    int movableTotal = 0;

    for (std::size_t iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      const auto& species = mTuneRuntime->species[iSpecies];

      if (!species.charged || species.etaCentralGivenActivity[info.mActivityClass].probability.empty()) {
        continue;
      }

      movableTotal += otherCounts[iSpecies];
    }

    if (movableTotal <= 0) {
      throw std::runtime_error("Ditto: cannot repair central charged coverage for the selected composition template");
    }

    std::uniform_int_distribution<int> moveDistribution(0, movableTotal - 1);

    int selected = moveDistribution(mRng);

    bool moved = false;

    for (std::size_t iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      const auto& species = mTuneRuntime->species[iSpecies];

      if (!species.charged || species.etaCentralGivenActivity[info.mActivityClass].probability.empty()) {
        continue;
      }

      const int count = otherCounts[iSpecies];

      if (selected < count) {
        --otherCounts[iSpecies];
        ++centralCounts[iSpecies];
        moved = true;
        break;
      }

      selected -= count;
    }

    if (!moved) {
      throw std::runtime_error("Ditto: failed to repair central charged composition");
    }

    --missingCentral;
  }

  int finalCentral = 0;
  int finalSelected = 0;

  for (std::size_t iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    finalCentral += centralCounts[iSpecies];

    finalSelected += centralCounts[iSpecies] + otherCounts[iSpecies];
  }

  if (finalCentral != info.mConditioningNch || finalSelected != info.mConditioningNSelected) {
    throw std::runtime_error("Ditto: repaired composition does not satisfy exact event multiplicities");
  }
}

double Generator::sampleAlias(const TuneAliasSampler& sampler)
{
  const std::size_t nBins = sampler.probability.size();

  if (nBins == 0) {
    throw std::runtime_error("Ditto: attempted to sample an empty alias table");
  }

  // Use one 64-bit RNG word for the complete histogram sample.
  //
  // The upper 32 bits define a discretized uniform variable used for both the
  // alias column and its Bernoulli decision. The lower 32 bits select the
  // continuous position inside the chosen histogram bin.
  const std::uint64_t random = mRng();

  const std::uint32_t selector = static_cast<std::uint32_t>(random >> 32);

  const std::uint64_t scaled = static_cast<std::uint64_t>(selector) * static_cast<std::uint64_t>(nBins);

  const std::size_t column = static_cast<std::size_t>(scaled >> 32);

  const double fraction = static_cast<double>(static_cast<std::uint32_t>(scaled)) * 0x1.0p-32;

  const std::size_t index = fraction < sampler.probability[column]
                              ? column
                              : static_cast<std::size_t>(sampler.alias[column]);

  const double withinBin = static_cast<double>(static_cast<std::uint32_t>(random)) * 0x1.0p-32;

  return sampler.lowEdges[index] + sampler.widths[index] * withinBin;
}

int Generator::sampleDiscreteAlias(const TuneDiscreteAliasSampler& sampler)
{
  const std::size_t nBins = sampler.probability.size();

  if (nBins == 0) {
    return -1;
  }

  const std::uint32_t selector = static_cast<std::uint32_t>(mRng() >> 32);

  const std::uint64_t scaled = static_cast<std::uint64_t>(selector) * static_cast<std::uint64_t>(nBins);

  const std::size_t column = static_cast<std::size_t>(scaled >> 32);

  const double fraction = static_cast<double>(static_cast<std::uint32_t>(scaled)) * 0x1.0p-32;

  const std::size_t index = fraction < sampler.probability[column]
                              ? column
                              : static_cast<std::size_t>(sampler.alias[column]);

  return sampler.values[index];
}

double Generator::samplePhi()
{
  return uniform(0.0, kTwoPi);
}

void Generator::makeParticle(TParticle* particle,
                             int pdg,
                             double mass,
                             double mass2,
                             double ptScale,
                             const TuneAliasSampler& ptSampler,
                             const TuneAliasSampler& etaSampler)
{
  using Clock = std::chrono::steady_clock;

  double pt = 0.0;
  double eta = 0.0;
  double phi = 0.0;

  if (detailedTimingEnabled()) {
    auto start = Clock::now();

    pt = sampleAlias(ptSampler) * ptScale;

    mTiming.mPtSampling += std::chrono::duration<double>(Clock::now() - start).count();

    start = Clock::now();

    eta = sampleAlias(etaSampler);

    mTiming.mEtaSampling += std::chrono::duration<double>(Clock::now() - start).count();

    start = Clock::now();

    phi = samplePhi();

    mTiming.mPhiSampling += std::chrono::duration<double>(Clock::now() - start).count();
  } else {
    pt = sampleAlias(ptSampler) * ptScale;
    eta = sampleAlias(etaSampler);
    phi = samplePhi();
  }

  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double e = 0.0;

  if (detailedTimingEnabled()) {
    const auto start = Clock::now();

    double sinPhi = 0.0;
    double cosPhi = 0.0;

#if defined(__GNUC__) || defined(__clang__)
    __builtin_sincos(phi,
                     &sinPhi,
                     &cosPhi);
#else
    sinPhi = std::sin(phi);
    cosPhi = std::cos(phi);
#endif

    px = pt * cosPhi;
    py = pt * sinPhi;
    pz = pt * std::sinh(eta);

    // px^2 + py^2 = pt^2 by construction.
    const double p2 = pt * pt + pz * pz;

    e = std::sqrt(p2 + mass2);

    mTiming.mMomentumMath += std::chrono::duration<double>(Clock::now() - start).count();
  } else {
    double sinPhi = 0.0;
    double cosPhi = 0.0;

#if defined(__GNUC__) || defined(__clang__)
    __builtin_sincos(phi,
                     &sinPhi,
                     &cosPhi);
#else
    sinPhi = std::sin(phi);
    cosPhi = std::cos(phi);
#endif

    px = pt * cosPhi;
    py = pt * sinPhi;
    pz = pt * std::sinh(eta);

    const double p2 = pt * pt + pz * pz;

    e = std::sqrt(p2 + mass2);
  }

  const auto setParticle = [&]() {
    particle->SetPdgCode(pdg);
    particle->SetStatusCode(mTune->mFinalStatus);

    particle->SetMother(0, -1);
    particle->SetMother(1, -1);

    particle->SetDaughter(0, -1);
    particle->SetDaughter(1, -1);

    particle->SetWeight(1.0f);
    particle->SetCalcMass(mass);

    particle->SetMomentum(px, py, pz, e);

    particle->SetProductionVertex(0.0, 0.0, 0.0, 0.0);

    particle->SetPolarisation(0.0, 0.0, 0.0);
  };

  if (detailedTimingEnabled()) {
    const auto start = Clock::now();

    setParticle();

    mTiming.mParticleSetters += std::chrono::duration<double>(Clock::now() - start).count();
  } else {
    setParticle();
  }
}

const EventInfo& Generator::generate()
{
  using Clock = std::chrono::steady_clock;

  const bool timing = timingEnabled();

  const bool detailed = detailedTimingEnabled();

  const auto generationStart = timing ? Clock::now() : Clock::time_point{};

  if (detailed) {
    auto start = Clock::now();

    mParticles->Clear("C");

    mTiming.mClearParticles += std::chrono::duration<double>(Clock::now() - start).count();

    start = Clock::now();

    mEvent = makeEventInfo();

    mTiming.mEventInfo += std::chrono::duration<double>(Clock::now() - start).count();
  } else {
    mParticles->Clear("C");
    mEvent = makeEventInfo();
  }

  int multiplicity = 0;

  if (detailed) {
    const auto start = Clock::now();

    multiplicity = sampleMultiplicity(mEvent);

    mTiming.mMultiplicity += std::chrono::duration<double>(Clock::now() - start).count();
  } else {
    multiplicity = sampleMultiplicity(mEvent);
  }

  if (multiplicity > mParticles->GetSize()) {
    if (detailed) {
      const auto start = Clock::now();

      mParticles->Expand(multiplicity);

      mTiming.mExpandArray += std::chrono::duration<double>(Clock::now() - start).count();
    } else {
      mParticles->Expand(multiplicity);
    }
  }

  const auto particleLoopStart = detailed ? Clock::now() : Clock::time_point{};

  int particleIndex = 0;

  const int nCentralCharged = mEvent.mConditioningNch;

  const int nOther = mEvent.mConditioningNSelected - mEvent.mConditioningNch;

  if (nCentralCharged < 0 || nOther < 0) {
    throw std::runtime_error("Ditto: invalid tuned event multiplicities");
  }

  if (detailed) {
    const auto start = Clock::now();

    sampleComposition(mEvent,
                      mTuneRuntime->currentCentralCounts,
                      mTuneRuntime->currentOtherCounts);

    mTiming.mSpeciesSampling += std::chrono::duration<double>(Clock::now() - start).count();
  } else {
    sampleComposition(mEvent,
                      mTuneRuntime->currentCentralCounts,
                      mTuneRuntime->currentOtherCounts);
  }

  const auto generateCounts = [&](const std::vector<int>& counts,
                                  bool centralCharged) {
    for (std::size_t iSpecies = 0; iSpecies < counts.size(); ++iSpecies) {
      const int count = counts[iSpecies];

      if (count <= 0) {
        continue;
      }

      const auto& species = mTuneRuntime->species[iSpecies];

      const auto& ptSampler = species.ptGivenActivity[mEvent.mActivityClass];

      const std::size_t nch = static_cast<std::size_t>(mEvent.mConditioningNch);

      const auto& ptScaleByNch = centralCharged ? species.centralPtScaleByNch
                                                : species.otherPtScaleByNch;

      if (nch >= ptScaleByNch.size()) {
        throw std::runtime_error("Ditto: exact-Nch pT correction index is out of range");
      }

      const double ptScale = ptScaleByNch[nch];

      const TuneAliasSampler* etaSampler = nullptr;

      if (centralCharged) {
        if (!species.charged) {
          throw std::runtime_error("Ditto: neutral species found in central-charged composition");
        }

        etaSampler = &species.etaCentralGivenActivity[mEvent.mActivityClass];
      } else if (species.charged) {
        etaSampler = &species.etaOutsideCentralGivenActivity[mEvent.mActivityClass];
      } else {
        etaSampler = &species.etaGivenActivity[mEvent.mActivityClass];
      }

      if (ptSampler.probability.empty() || !etaSampler || etaSampler->probability.empty()) {
        throw std::runtime_error("Ditto: empty tuned kinematic sampler for PDG " + std::to_string(species.pdg));
      }

      for (int i = 0; i < count;
           ++i) {
        TParticle* particle = nullptr;

        if (detailed) {
          const auto start = Clock::now();

          particle = static_cast<TParticle*>(mParticles->ConstructedAt(particleIndex));

          mTiming.mConstructedAt += std::chrono::duration<double>(Clock::now() - start).count();
        } else {
          particle = static_cast<TParticle*>(mParticles->ConstructedAt(particleIndex));
        }

        makeParticle(particle,
                     species.pdg,
                     species.mass,
                     species.mass2,
                     ptScale,
                     ptSampler,
                     *etaSampler);

        ++particleIndex;
      }
    }
  };

  // The empirical template already contains correlated species counts.
  // The central-coverage repair, when needed, preserves all species totals.
  generateCounts(mTuneRuntime->currentCentralCounts, true);

  generateCounts(mTuneRuntime->currentOtherCounts, false);

  if (particleIndex != multiplicity) {
    throw std::runtime_error("Ditto: constrained species counts do not reproduce Nselected");
  }

  // Overlay the signal only after the minimum-bias event is fully fixed. The
  // injector consumes the same activity class selected from the sampled exact
  // Nch; it never resamples activity and does not alter the conditioning
  // variables of the background event.
  mEvent.mInjectedParticles = 0;
  mEvent.mInjectedDecayProducts = 0;
  if (mSignalInjector) {
    const auto signal = mSignalInjector->inject(*mParticles,
                                                particleIndex,
                                                mEvent.mActivityClass);
    mEvent.mInjectedParticles = signal.mInjectedParticles;
    mEvent.mInjectedDecayProducts = signal.mDecayProducts;
    particleIndex += signal.totalParticles();
  }

  if (detailed) {
    mTiming.mParticleLoop += std::chrono::duration<double>(Clock::now() - particleLoopStart).count();
  }

  // Stop pure-generation timing before ROOT I/O.
  if (timing) {
    const auto generationStop = Clock::now();

    mTiming.mGeneration += std::chrono::duration<double>(generationStop - generationStart).count();

    ++mTiming.mGeneratedEvents;

    mTiming.mGeneratedParticles += static_cast<std::uint64_t>(particleIndex);
  }

  if (mTree) {
    if (timing) {
      const auto start = Clock::now();

      mTree->Fill();

      mTiming.mTreeFill += std::chrono::duration<double>(Clock::now() - start).count();

      ++mTiming.mTreeFills;
    } else {
      mTree->Fill();
    }
  }

  ++mEventCounter;

  return mEvent;
}

void Generator::loadParticles(Pythia8::Event& event, bool reset) const
{
  using Clock = std::chrono::steady_clock;
  const bool timing = timingEnabled();
  const bool detailed = detailedTimingEnabled();
  const auto exportStart = timing ? Clock::now() : Clock::time_point{};

  if (reset) {
    if (detailed) {
      const auto start = Clock::now();
      event.reset();
      mTiming.mPythiaReset += std::chrono::duration<double>(Clock::now() - start).count();
    } else {
      event.reset();
    }
  }

  const int baseIndex = event.size();
  const int nParticles = mParticles->GetEntriesFast();

  const auto appendParticles = [&]() {
    for (int i = 0; i < nParticles; ++i) {
      const auto* p = static_cast<const TParticle*>(mParticles->UncheckedAt(i));
      if (!p) {
        continue;
      }

      const int m1 = p->GetMother(0);
      const int m2 = p->GetMother(1);
      const int d1 = p->GetDaughter(0);
      const int d2 = p->GetDaughter(1);

      const int mother1 = m1 < 0 ? 0 : baseIndex + m1;
      const int mother2 = m2 < 0 ? 0 : baseIndex + m2;
      const int daughter1 = d1 < 0 ? 0 : baseIndex + d1;
      const int daughter2 = d2 < 0 ? 0 : baseIndex + d2;

      const int index = event.append(p->GetPdgCode(),
                                     p->GetStatusCode(),
                                     mother1,
                                     mother2,
                                     daughter1,
                                     daughter2,
                                     0,
                                     0,
                                     p->Px(),
                                     p->Py(),
                                     p->Pz(),
                                     p->Energy(),
                                     p->GetCalcMass());

      event[index].vProd(p->Vx(), p->Vy(), p->Vz(), p->T());
    }
  };

  if (detailed) {
    const auto start = Clock::now();
    appendParticles();
    mTiming.mPythiaAppend += std::chrono::duration<double>(Clock::now() - start).count();
  } else {
    appendParticles();
  }

  // When we created a fresh PYTHIA event, update entry 0 to represent the sum
  // of the Ditto particles. If appending to an existing event, leave the
  // existing event-level pseudo-particle untouched.
  if (reset && event.size() > 0) {
    const auto updateSystem = [&]() {
      double px = 0.0;
      double py = 0.0;
      double pz = 0.0;
      double e = 0.0;

      for (int i = 0; i < nParticles; ++i) {
        const auto* p = static_cast<const TParticle*>(mParticles->UncheckedAt(i));
        if (!p) {
          continue;
        }

        px += p->Px();
        py += p->Py();
        pz += p->Pz();
        e += p->Energy();
      }

      const double m2 = e * e - px * px - py * py - pz * pz;
      event[0].p(px, py, pz, e);
      event[0].m(std::sqrt(std::max(0.0, m2)));
    };

    if (detailed) {
      const auto start = Clock::now();
      updateSystem();
      mTiming.mPythiaSystemSum += std::chrono::duration<double>(Clock::now() - start).count();
    } else {
      updateSystem();
    }
  }

  if (timing) {
    mTiming.mPythiaExport += std::chrono::duration<double>(Clock::now() - exportStart).count();
    ++mTiming.mPythiaExports;
  }
}

void Generator::registerTTreeOutput(const std::string& fileName)
{
  if (mOutputFile || mTree) {
    throw std::logic_error("Ditto: TTree output is already registered");
  }

  mOutputFile = TFile::Open(fileName.c_str(), "RECREATE");
  if (!mOutputFile || mOutputFile->IsZombie()) {
    delete mOutputFile;
    mOutputFile = nullptr;
    throw std::runtime_error("Ditto: could not create ROOT file " + fileName);
  }

  mTree = new TTree(std::string(kTreeName).c_str(), "Ditto generated events");
  mTree->Branch(std::string(kBranchName).c_str(), &mParticles, 256000, 0);
  mTree->SetDirectory(mOutputFile);
}

} // namespace Ditto
