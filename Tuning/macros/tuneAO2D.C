#include "DittoAO2DTuner.h"

#include <fstream>

void tuneAO2D(const std::string& inputFile = "aodlist.txt",
              const std::string& refPythiaCard = "Tuning/cards/pythia8_inel_136tev.cfg")
{
  Ditto::AO2DTunerConfig cfg;

  // Read the input file list.
  std::ifstream file(inputFile);
  std::string line;
  while (std::getline(file, line)) {
    if (!line.empty()) {
      cfg.mInputFiles.push_back(line);
    }
  }

  cfg.mPythiaCard = refPythiaCard;

  // 0 = all MC collisions in all input files.
  cfg.mMaxEvents = 1E9;
  cfg.mProgressEvery = 10000;

  // Closest AO2D equivalent to PYTHIA isFinal().
  cfg.mParticleSelection = Ditto::AO2DParticleSelection::GeneratorFinal;

  // The normal TuneAccumulator configuration is available here too.
  cfg.mActivityEtaMax = 0.5;
  cfg.mParticleEtaMax = 5.0;
  cfg.mPtMax = 20.0;

  Ditto::AO2DTuner tuner(cfg);
  tuner.run();
  tuner.save("Ditto_tune_AO2D.root");
}
