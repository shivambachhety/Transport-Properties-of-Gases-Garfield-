#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <map>

#include <TMath.h>

#include "MediumMagboltz.hh"
#include "MagboltzInterface.hh"
#include "Random.hh"
#include "FundamentalConstants.hh"
#include "GarfieldConstants.hh"
#include "OpticalData.hh"

namespace Garfield {

const int MediumMagboltz::DxcTypeRad = 0;
const int MediumMagboltz::DxcTypeCollIon = 1;
const int MediumMagboltz::DxcTypeCollNonIon = -1;

MediumMagboltz::MediumMagboltz()
    : MediumGas(),
      m_eFinal(40.),
      m_eStep(m_eFinal / nEnergySteps),
      m_eHigh(1.e4),
      m_eHighLog(log(m_eHigh)),
      m_lnStep(1.),
      m_useAutoAdjust(true),
      m_useCsOutput(false),
      m_nTerms(0),
      m_useAnisotropic(true),
      m_nPenning(0),
      m_useDeexcitation(false),
      m_useRadTrap(true),
      m_useOpalBeaty(true),
      m_useGreenSawada(false),
      m_eFinalGamma(20.),
      m_eStepGamma(m_eFinalGamma / nEnergyStepsGamma) {
 
  fit3d4p = fitHigh4p = 1.;
  fit3dQCO2 = fit3dQCH4 = fit3dQC2H6 = 1.;
  fit3dEtaCO2 = fit3dEtaCH4 = fit3dEtaC2H6 = 0.5;
  fit4pEtaCH4 = fit4pEtaC2H6 = 0.5;
  fit4sEtaC2H6 = 0.5;
  fitLineCut = 1000;

  m_className = "MediumMagboltz";

  // Set physical constants in Magboltz common blocks.
  Magboltz::cnsts_.echarg = ElementaryCharge * 1.e-15;
  Magboltz::cnsts_.emass = ElectronMassGramme;
  Magboltz::cnsts_.amu = AtomicMassUnit;
  Magboltz::cnsts_.pir2 = BohrRadius * BohrRadius * Pi;
  Magboltz::inpt_.ary = RydbergEnergy;

  // Set parameters in Magboltz common blocks.
  Magboltz::inpt_.nGas = m_nComponents;
  Magboltz::inpt_.nStep = nEnergySteps;
  // Select the scattering model.
  Magboltz::inpt_.nAniso = 2;
  // Max. energy [eV]
  Magboltz::inpt_.efinal = m_eFinal;
  // Energy step size [eV]
  Magboltz::inpt_.estep = m_eStep;
  // Temperature and pressure
  Magboltz::inpt_.akt = BoltzmannConstant * m_temperature;
  Magboltz::inpt_.tempc = m_temperature - ZeroCelsius;
  Magboltz::inpt_.torr = m_pressure;
  // Disable Penning transfer.
  Magboltz::inpt_.ipen = 0;

  // Initialise Penning parameters
  for (int i = nMaxLevels; i--;) {
    m_rPenning[i] = 0.;
    m_lambdaPenning[i] = 0.;
  }

  m_isChanged = true;

  EnableDrift();
  EnablePrimaryIonisation();
  m_microscopic = true;

  // Initialize the collision counters.
  m_nCollisionsDetailed.clear();
  for (int i = nCsTypes; i--;) m_nCollisions[i] = 0;
  for (int i = nCsTypesGamma; i--;) m_nPhotonCollisions[i] = 0;

  m_ionProducts.clear();
  m_dxcProducts.clear();

  for (unsigned int i = 0; i < m_nMaxGases; ++i) m_scaleExc[i] = 1.;
}

bool MediumMagboltz::SetMaxElectronEnergy(const double e) {

  if (e <= Small) {
    std::cerr << m_className << "::SetMaxElectronEnergy:\n";
    std::cerr << "    Provided upper electron energy limit (" << e
              << " eV) is too small.\n";
    return false;
  }
  m_eFinal = e;

  // Determine the energy interval size.
  if (m_eFinal <= m_eHigh) {
    m_eStep = m_eFinal / nEnergySteps;
  } else {
    m_eStep = m_eHigh / nEnergySteps;
  }

  // Set max. energy and step size also in Magboltz common block.
  Magboltz::inpt_.efinal = m_eFinal;
  Magboltz::inpt_.estep = m_eStep;

  // Force recalculation of the scattering rates table.
  m_isChanged = true;

  return true;
}

bool MediumMagboltz::SetMaxPhotonEnergy(const double e) {

  if (e <= Small) {
    std::cerr << m_className << "::SetMaxPhotonEnergy:\n";
    std::cerr << "    Provided upper photon energy limit (" << e
              << " eV) is too small.\n";
    return false;
  }
  m_eFinalGamma = e;

  // Determine the energy interval size.
  m_eStepGamma = m_eFinalGamma / nEnergyStepsGamma;

  // Force recalculation of the scattering rates table.
  m_isChanged = true;

  return true;
}

void MediumMagboltz::SetSplittingFunctionOpalBeaty() {

  m_useOpalBeaty = true;
  m_useGreenSawada = false;
}

void MediumMagboltz::SetSplittingFunctionGreenSawada() {

  m_useOpalBeaty = false;
  m_useGreenSawada = true;
  if (m_isChanged) return;

  bool allset = true;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (!m_hasGreenSawada[i]) {
      if (allset) {
        std::cout << m_className << "::SetSplittingFunctionGreenSawada:\n";
        allset = false;
      }
      std::cout << "    Fit parameters for " << m_gas[i] << " not available.\n";
      std::cout << "    Opal-Beaty formula is used instead.\n";
    }
  }
}

void MediumMagboltz::SetSplittingFunctionFlat() {

  m_useOpalBeaty = false;
  m_useGreenSawada = false;
}

void MediumMagboltz::EnableDeexcitation() {

  if (m_usePenning) {
    std::cout << m_className << "::EnableDeexcitation:\n";
    std::cout << "    Penning transfer will be switched off.\n";
  }
  // if (m_useRadTrap) {
  //   std::cout << "    Radiation trapping is switched on.\n";
  // } else {
  //   std::cout << "    Radiation trapping is switched off.\n";
  // }
  m_usePenning = false;
  m_useDeexcitation = true;
  m_isChanged = true;
  m_dxcProducts.clear();
}

void MediumMagboltz::EnableRadiationTrapping() {

  m_useRadTrap = true;
  if (!m_useDeexcitation) {
    std::cout << m_className << "::EnableRadiationTrapping:\n";
    std::cout << "    Radiation trapping is enabled"
              << " but de-excitation is not.\n";
  } else {
    m_isChanged = true;
  }
}

void MediumMagboltz::EnablePenningTransfer(const double r,
                                           const double lambda) {

  if (r < 0. || r > 1.) {
    std::cerr << m_className << "::EnablePenningTransfer:\n";
    std::cerr << "    Penning transfer probability must be "
              << " in the range [0, 1].\n";
    return;
  }

  m_rPenningGlobal = r;
  if (lambda < Small) {
    m_lambdaPenningGlobal = 0.;
  } else {
    m_lambdaPenningGlobal = lambda;
  }

  std::cout << m_className << "::EnablePenningTransfer:\n";
  std::cout << "    Global Penning transfer parameters set to: \n";
  std::cout << "    r      = " << m_rPenningGlobal << "\n";
  std::cout << "    lambda = " << m_lambdaPenningGlobal << " cm\n";

  for (unsigned int i = 0; i < m_nTerms; ++i) {
    m_rPenning[i] = m_rPenningGlobal;
    m_lambdaPenning[i] = m_lambdaPenningGlobal;
  }

  if (m_useDeexcitation) {
    std::cout << m_className << "::EnablePenningTransfer:\n";
    std::cout << "    Deexcitation handling will be switched off.\n";
  }
  m_usePenning = true;
}

void MediumMagboltz::EnablePenningTransfer(const double r, const double lambda,
                                           std::string gasname) {

  if (r < 0. || r > 1.) {
    std::cerr << m_className << "::EnablePenningTransfer:\n";
    std::cerr << "    Penning transfer probability must be "
              << " in the range [0, 1].\n";
    return;
  }

  // Get the "standard" name of this gas.
  if (!GetGasName(gasname, gasname)) {
    std::cerr << m_className << "::EnablePenningTransfer:\n";
    std::cerr << "    Unknown gas name.\n";
    return;
  }

  // Look for this gas in the present gas mixture.
  bool found = false;
  int iGas = -1;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == gasname) {
      m_rPenningGas[i] = r;
      if (lambda < Small) {
        m_lambdaPenningGas[i] = 0.;
      } else {
        m_lambdaPenningGas[i] = lambda;
      }
      found = true;
      iGas = i;
      break;
    }
  }

  if (!found) {
    std::cerr << m_className << "::EnablePenningTransfer:\n";
    std::cerr << "    Specified gas (" << gasname
              << ") is not part of the present gas mixture.\n";
    return;
  }

  // Make sure that the collision rate table is updated.
  if (m_isChanged) {
    if (!Mixer()) {
      std::cerr << m_className << "::EnablePenningTransfer:\n";
      std::cerr << "    Error calculating the collision rates table.\n";
      return;
    }
    m_isChanged = false;
  }

  unsigned int nLevelsFound = 0;
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    if (int(m_csType[i] / nCsTypes) != iGas) continue;
    if (m_csType[i] % nCsTypes == ElectronCollisionTypeExcitation) {
      ++nLevelsFound;
    }
    m_rPenning[i] = m_rPenningGas[iGas];
    m_lambdaPenning[i] = m_lambdaPenningGas[iGas];
  }

  if (nLevelsFound > 0) {
    std::cout << m_className << "::EnablePenningTransfer:\n";
    std::cout << "    Penning transfer parameters for " << nLevelsFound
              << " excitation levels set to:\n";
    std::cout << "      r      = " << m_rPenningGas[iGas] << "\n";
    std::cout << "      lambda = " << m_lambdaPenningGas[iGas] << " cm\n";
  } else {
    std::cerr << m_className << "::EnablePenningTransfer:\n";
    std::cerr << "    Specified gas (" << gasname
              << ") has no excitation levels in the present energy range.\n";
  }

  m_usePenning = true;
}

void MediumMagboltz::DisablePenningTransfer() {

  for (unsigned int i = 0; i < m_nTerms; ++i) {
    m_rPenning[i] = 0.;
    m_lambdaPenning[i] = 0.;
  }
  m_rPenningGlobal = 0.;
  m_lambdaPenningGlobal = 0.;

  for (unsigned int i = 0; i < m_nMaxGases; ++i) {
    m_rPenningGas[i] = 0.;
    m_lambdaPenningGas[i] = 0.;
  }

  m_usePenning = false;
}

void MediumMagboltz::DisablePenningTransfer(std::string gasname) {

  // Get the "standard" name of this gas.
  if (!GetGasName(gasname, gasname)) {
    std::cerr << m_className << "::DisablePenningTransfer:\n";
    std::cerr << "    Gas " << gasname << " is not defined.\n";
    return;
  }

  // Look for this gas in the present gas mixture.
  bool found = false;
  int iGas = -1;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == gasname) {
      m_rPenningGas[i] = 0.;
      m_lambdaPenningGas[i] = 0.;
      found = true;
      iGas = i;
      break;
    }
  }

  if (!found) {
    std::cerr << m_className << "::DisablePenningTransfer:\n";
    std::cerr << "    Specified gas (" << gasname
              << ") is not part of the present gas mixture.\n";
    return;
  }

  unsigned int nLevelsFound = 0;
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    if (int(m_csType[i] / nCsTypes) == iGas) {
      m_rPenning[i] = 0.;
      m_lambdaPenning[i] = 0.;
    } else {
      if (m_csType[i] % nCsTypes == ElectronCollisionTypeExcitation &&
          m_rPenning[i] > Small) {
        ++nLevelsFound;
      }
    }
  }

  if (nLevelsFound == 0) {
    // There are no more excitation levels with r > 0.
    std::cout << m_className << "::DisablePenningTransfer:\n"
              << "    Penning transfer globally switched off.\n";
    m_usePenning = false;
  }
}

void MediumMagboltz::SetExcitationScalingFactor(const double r,
                                                std::string gasname) {

  if (r <= 0.) {
    std::cerr << m_className << "::SetScalingFactor:\n";
    std::cerr << "    Incorrect value for scaling factor: " << r << "\n";
    return;
  }

  // Get the "standard" name of this gas.
  if (!GetGasName(gasname, gasname)) {
    std::cerr << m_className << "::SetExcitationScalingFactor:\n";
    std::cerr << "    Unknown gas name.\n";
    return;
  }

  // Look for this gas in the present gas mixture.
  bool found = false;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == gasname) {
      m_scaleExc[i] = r;
      found = true;
      break;
    }
  }

  if (!found) {
    std::cerr << m_className << "::SetExcitationScalingFactor:\n";
    std::cerr << "    Specified gas (" << gasname
              << ") is not part of the present gas mixture.\n";
    return;
  }

  // Make sure that the collision rate table is updated.
  m_isChanged = true;
}

bool MediumMagboltz::Initialise(const bool verbose) {

  if (!m_isChanged) {
    if (m_debug) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Nothing changed.\n";
    }
    return true;
  }
  if (!Mixer(verbose)) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Error calculating the collision rates table.\n";
    return false;
  }
  m_isChanged = false;
  return true;
}

void MediumMagboltz::PrintGas() {

  MediumGas::PrintGas();

  if (m_isChanged) {
    if (!Initialise()) return;
  }

  std::cout << m_className << "::PrintGas:\n";
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    // Collision type
    int type = m_csType[i] % nCsTypes;
    int ngas = int(m_csType[i] / nCsTypes);
    // Description (from Magboltz)
    std::string descr = std::string(50, ' ');
    for (int j = 50; j--;) descr[j] = m_description[i][j];
    // Threshold energy
    double e = m_rgas[ngas] * m_energyLoss[i];
    std::cout << "    Level " << i << ": " << descr << "\n";
    std::cout << "        Type " << type;
    if (type == ElectronCollisionTypeElastic) {
      std::cout << " (elastic)\n";
    } else if (type == ElectronCollisionTypeIonisation) {
      std::cout << " (ionisation)\n";
      std::cout << "        Ionisation threshold: " << e << " eV\n";
    } else if (type == ElectronCollisionTypeAttachment) {
      std::cout << " (attachment)\n";
    } else if (type == ElectronCollisionTypeInelastic) {
      std::cout << " (inelastic)\n";
      std::cout << "        Energy loss: " << e << " eV\n";
    } else if (type == ElectronCollisionTypeExcitation) {
      std::cout << " (excitation)\n";
      std::cout << "        Excitation energy: " << e << " eV\n";
    } else if (type == ElectronCollisionTypeSuperelastic) {
      std::cout << " (super-elastic)\n";
      std::cout << "        Energy gain: " << -e << " eV\n";
    } else {
      std::cout << " (unknown)\n";
    }
    if (type == ElectronCollisionTypeExcitation && m_usePenning &&
        e > m_minIonPot) {
      std::cout << "        Penning transfer coefficient: " << m_rPenning[i]
                << "\n";
    } else if (type == ElectronCollisionTypeExcitation && m_useDeexcitation) {
      const int idxc = m_iDeexcitation[i];
      if (idxc < 0 || idxc >= (int)m_deexcitations.size()) {
        std::cout << "        Deexcitation cascade not implemented.\n";
        continue;
      }
      if (m_deexcitations[idxc].osc > 0.) {
        std::cout << "        Oscillator strength: " << m_deexcitations[idxc].osc
                  << "\n";
      }
      std::cout << "        Decay channels:\n";
      for (int j = 0; j < m_deexcitations[idxc].nChannels; ++j) {
        if (m_deexcitations[idxc].type[j] == DxcTypeRad) {
          std::cout << "          Radiative decay to ";
          if (m_deexcitations[idxc].final[j] < 0) {
            std::cout << "ground state: ";
          } else {
            std::cout << m_deexcitations[m_deexcitations[idxc].final[j]].label
                      << ": ";
          }
        } else if (m_deexcitations[idxc].type[j] == DxcTypeCollIon) {
          if (m_deexcitations[idxc].final[j] < 0) {
            std::cout << "          Penning ionisation: ";
          } else {
            std::cout << "          Associative ionisation: ";
          }
        } else if (m_deexcitations[idxc].type[j] == DxcTypeCollNonIon) {
          if (m_deexcitations[idxc].final[j] >= 0) {
            std::cout << "          Collision-induced transition to "
                      << m_deexcitations[m_deexcitations[idxc].final[j]].label
                      << ": ";
          } else {
            std::cout << "          Loss: ";
          }
        }
        if (j == 0) {
          std::cout << std::setprecision(5) << m_deexcitations[idxc].p[j] * 100.
                    << "%\n";
        } else {
          std::cout << std::setprecision(5) << (m_deexcitations[idxc].p[j] -
                                                m_deexcitations[idxc].p[j - 1]) *
                                                   100. << "%\n";
        }
      }
    }
  }
}

double MediumMagboltz::GetElectronNullCollisionRate(const int band) {

  // If necessary, update the collision rates table.
  if (m_isChanged) {
    if (!Mixer()) {
      std::cerr << m_className << "::GetElectronNullCollisionRate:\n";
      std::cerr << "     Error calculating the collision rates table.\n";
      return 0.;
    }
    m_isChanged = false;
  }

  if (m_debug && band > 0) {
    std::cerr << m_className << "::GetElectronNullCollisionRate:\n";
    std::cerr << "    Warning: unexpected band index.\n";
  }

  return m_cfNull;
}

double MediumMagboltz::GetElectronCollisionRate(const double e,
                                                const int band) {

  // Check if the electron energy is within the currently set range.
  if (e <= 0.) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n";
    std::cerr << "    Electron energy must be greater than zero.\n";
    return m_cfTot[0];
  }
  if (e > m_eFinal && m_useAutoAdjust) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n";
    std::cerr << "    Collision rate at " << e
              << " eV is not included in the current table.\n";
    std::cerr << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  }

  // If necessary, update the collision rates table.
  if (m_isChanged) {
    if (!Mixer()) {
      std::cerr << m_className << "::GetElectronCollisionRate:\n";
      std::cerr << "    Error calculating the collision rates table.\n";
      return 0.;
    }
    m_isChanged = false;
  }

  if (m_debug && band > 0) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n";
    std::cerr << "    Warning: unexpected band index.\n";
  }

  // Get the energy interval.
  int iE = 0;
  if (e <= m_eHigh) {
    // Linear binning
    iE = int(e / m_eStep);
    if (iE >= nEnergySteps) return m_cfTot[nEnergySteps - 1];
    if (iE < 0) return m_cfTot[0];
    return m_cfTot[iE];
  }

  // Logarithmic binning
  const double eLog = log(e);
  iE = int((eLog - m_eHighLog) / m_lnStep);
  // Calculate the collision rate by log-log interpolation.
  const double fmax = m_cfTotLog[iE];
  const double fmin = iE == 0 ? log(m_cfTot[nEnergySteps - 1]) : m_cfTotLog[iE - 1];
  const double emin = m_eHighLog + iE * m_lnStep;
  const double f = fmin + (eLog - emin) * (fmax - fmin) / m_lnStep;
  return exp(f);
}

double MediumMagboltz::GetElectronCollisionRate(const double e, 
                                                const unsigned int level,
                                                const int band) {

  // Check if the electron energy is within the currently set range.
  if (e <= 0.) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n";
    std::cerr << "    Electron energy must be greater than zero.\n";
    return 0.;
  }

  // Check if the level exists.
  if (level >= m_nTerms) {
    std::cerr << m_className << "::GetElectronCollisionRate:\n";
    std::cerr << "    Level " << level << " does not exist.\n";
    std::cerr << "    The present gas mixture has " << m_nTerms
              << " cross-section terms.\n";
    return 0.;
  }

  // Get the total scattering rate.
  double rate = GetElectronCollisionRate(e, band);
  // Get the energy interval.
  int iE = 0;
  if (e <= m_eHigh) {
    // Linear binning
    iE = int(e / m_eStep);
    if (iE >= nEnergySteps) return m_cfTot[nEnergySteps - 1];
    if (level == 0) {
      rate *= m_cf[iE][0];
    } else {
      rate *= m_cf[iE][level] - m_cf[iE][level - 1];
    }
  } else {
    // Logarithmic binning
    iE = int((log(e) - m_eHighLog) / m_lnStep);
    if (level == 0) {
      rate *= m_cfLog[iE][0];
    } else {
      rate *= m_cfLog[iE][level] - m_cfLog[iE][level - 1];
    }
  }
  return rate;
}

bool MediumMagboltz::GetElectronCollision(const double e, int& type, int& level,
                                          double& e1, double& dx, double& dy,
                                          double& dz, int& nion, int& ndxc,
                                          int& band) {

  // Check if the electron energy is within the currently set range.
  if (e > m_eFinal && m_useAutoAdjust) {
    std::cerr << m_className << "::GetElectronCollision:\n";
    std::cerr << "    Provided electron energy  (" << e
              << " eV) exceeds current energy range  (" << m_eFinal << " eV).\n";
    std::cerr << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  } else if (e <= 0.) {
    std::cerr << m_className << "::GetElectronCollision:\n";
    std::cerr << "    Electron energy must be greater than zero.\n";
    return false;
  }

  // If necessary, update the collision rates table.
  if (m_isChanged) {
    if (!Mixer()) {
      std::cerr << m_className << "::GetElectronCollision:\n";
      std::cerr << "    Error calculating the collision rates table.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (m_debug && band > 0) {
    std::cerr << m_className << "::GetElectronCollision:\n";
    std::cerr << "    Warning: unexpected band index.\n";
  }

  double angCut = 1.;
  double angPar = 0.5;

  if (e <= m_eHigh) {
    // Linear binning
    // Get the energy interval.
    int iE = int(e / m_eStep);
    if (iE >= nEnergySteps) iE = nEnergySteps - 1;
    if (iE < 0) iE = 0;

    // Sample the scattering process.
    const double r = RndmUniform();
    int iLow = 0;
    int iUp = m_nTerms - 1;
    if (r <= m_cf[iE][iLow]) {
      level = iLow;
    } else if (r >= m_cf[iE][iUp]) {
      level = iUp;
    } else {
      int iMid;
      while (iUp - iLow > 1) {
        iMid = (iLow + iUp) >> 1;
        if (r < m_cf[iE][iMid]) {
          iUp = iMid;
        } else {
          iLow = iMid;
        }
      }
      level = iUp;
    }
    // Get the angular distribution parameters.
    angCut = m_scatCut[iE][level];
    angPar = m_scatParameter[iE][level];
  } else {
    // Logarithmic binning
    // Get the energy interval.
    int iE = int(log(e / m_eHigh) / m_lnStep);
    if (iE < 0) iE = 0;
    if (iE >= nEnergyStepsLog) iE = nEnergyStepsLog - 1;
    // Sample the scattering process.
    const double r = RndmUniform();
    int iLow = 0;
    int iUp = m_nTerms - 1;
    if (r <= m_cfLog[iE][iLow]) {
      level = iLow;
    } else if (r >= m_cfLog[iE][iUp]) {
      level = iUp;
    } else {
      int iMid;
      while (iUp - iLow > 1) {
        iMid = (iLow + iUp) >> 1;
        if (r < m_cfLog[iE][iMid]) {
          iUp = iMid;
        } else {
          iLow = iMid;
        }
      }
      level = iUp;
    }
    // Get the angular distribution parameters.
    angCut = m_scatCutLog[iE][level];
    angPar = m_scatParameterLog[iE][level];
  }

  // Extract the collision type.
  type = m_csType[level] % nCsTypes;
  const int igas = int(m_csType[level] / nCsTypes);
  // Increase the collision counters.
  ++m_nCollisions[type];
  ++m_nCollisionsDetailed[level];

  // Get the energy loss for this process.
  double loss = m_energyLoss[level];
  nion = ndxc = 0;

  if (type == ElectronCollisionTypeIonisation) {
    // Sample the secondary electron energy according to
    // the Opal-Beaty-Peterson parameterisation.
    double esec = 0.;
    if (m_useOpalBeaty) {
      // Get the splitting parameter.
      const double w = m_wOpalBeaty[level];
      esec = w * tan(RndmUniform() * atan(0.5 * (e - loss) / w));
      // Rescaling (SST)
      // esec = w * pow(esec / w, 0.9524);
    } else if (m_useGreenSawada) {
      const double w = m_gsGreenSawada[igas] * e / (e + m_gbGreenSawada[igas]);
      const double esec0 =
          m_tsGreenSawada[igas] - m_taGreenSawada[igas] / (e + m_tbGreenSawada[igas]);
      const double r = RndmUniform();
      esec = esec0 + w * tan((r - 1.) * atan(esec0 / w) +
                             r * atan((0.5 * (e - loss) - esec0) / w));
    } else {
      esec = RndmUniform() * (e - loss);
    }
    if (esec <= 0) esec = Small;
    loss += esec;
    m_ionProducts.clear();
    // Add the secondary electron.
    ionProd newIonProd;
    newIonProd.type = IonProdTypeElectron;
    newIonProd.energy = esec;
    m_ionProducts.push_back(newIonProd);
    // Add the ion.
    newIonProd.type = IonProdTypeIon;
    newIonProd.energy = 0.;
    m_ionProducts.push_back(newIonProd);
    nion = 2;
  } else if (type == ElectronCollisionTypeExcitation) {
    // if (m_gas[igas] == "CH4" && loss * m_rgas[igas] < 13.35 && e > 12.65) {
    //   if (RndmUniform() < 0.5) {
    //     loss = 8.55 + RndmUniform() * (13.3 - 8.55);
    //     loss /= m_rgas[igas];
    //   } else {
    //     loss = std::max(Small, RndmGaussian(loss * m_rgas[igas], 1.));
    //     loss /= m_rgas[igas];
    //   }
    // }
    // Follow the de-excitation cascade (if switched on).
    if (m_useDeexcitation && m_iDeexcitation[level] >= 0) {
      int fLevel = 0;
      ComputeDeexcitationInternal(m_iDeexcitation[level], fLevel);
      ndxc = m_dxcProducts.size();
    } else if (m_usePenning) {
      m_dxcProducts.clear();
      // Simplified treatment of Penning ionisation.
      // If the energy threshold of this level exceeds the
      // ionisation potential of one of the gases,
      // create a new electron (with probability m_rPenning).
      if (m_energyLoss[level] * m_rgas[igas] > m_minIonPot &&
          RndmUniform() < m_rPenning[level]) {
        // The energy of the secondary electron is assumed to be given by
        // the difference of excitation and ionisation threshold.
        double esec = m_energyLoss[level] * m_rgas[igas] - m_minIonPot;
        if (esec <= 0) esec = Small;
        // Add the secondary electron to the list.
        dxcProd newDxcProd;
        newDxcProd.t = 0.;
        newDxcProd.s = 0.;
        if (m_lambdaPenning[level] > Small) {
          // Uniform distribution within a sphere of radius lambda
          newDxcProd.s = m_lambdaPenning[level] * pow(RndmUniformPos(), 1. / 3.);
        }
        newDxcProd.energy = esec;
        newDxcProd.type = DxcProdTypeElectron;
        m_dxcProducts.push_back(newDxcProd);
        ndxc = 1;
        ++m_nPenning;
      }
    }
  }

  // Make sure the energy loss is smaller than the energy.
  if (e < loss) loss = e - 0.0001;

  // Determine the scattering angle.
  double ctheta0 = 1. - 2. * RndmUniform();
  if (m_useAnisotropic) {
    switch (m_scatModel[level]) {
      case 0:
        break;
      case 1:
        ctheta0 = 1. - RndmUniform() * angCut;
        if (RndmUniform() > angPar) ctheta0 = -ctheta0;
        break;
      case 2:
        ctheta0 = (ctheta0 + angPar) / (1. + angPar * ctheta0);
        break;
      default:
        std::cerr << m_className << "::GetElectronCollision:\n";
        std::cerr << "    Unknown scattering model. \n";
        std::cerr << "    Using isotropic distribution.\n";
        break;
    }
  }

  const double s1 = m_rgas[igas];
  const double s2 = (s1 * s1) / (s1 - 1.);
  const double theta0 = acos(ctheta0);
  const double arg = std::max(1. - s1 * loss / e, Small);
  const double d = 1. - ctheta0 * sqrt(arg);

  // Update the energy.
  e1 = std::max(e * (1. - loss / (s1 * e) - 2. * d / s2), Small);
  double q = std::min(sqrt((e / e1) * arg) / s1, 1.);
  const double theta = asin(q * sin(theta0));
  double ctheta = cos(theta);
  if (ctheta0 < 0.) {
    const double u = (s1 - 1.) * (s1 - 1.) / arg;
    if (ctheta0 * ctheta0 > u) ctheta = -ctheta;
  }
  const double stheta = sin(theta);
  // Calculate the direction after the collision.
  dz = std::min(dz, 1.);
  const double argZ = sqrt(dx * dx + dy * dy);

  // Azimuth is chosen at random.
  const double phi = TwoPi * RndmUniform();
  const double cphi = cos(phi);
  const double sphi = sin(phi);

  if (argZ == 0.) {
    dz = ctheta;
    dx = cphi * stheta;
    dy = sphi * stheta;
  } else {
    const double a = stheta / argZ;
    const double dz1 = dz * ctheta + argZ * stheta * sphi;
    const double dy1 = dy * ctheta + a * (dx * cphi - dy * dz * sphi);
    const double dx1 = dx * ctheta - a * (dy * cphi + dx * dz * sphi);
    dz = dz1;
    dy = dy1;
    dx = dx1;
  }

  return true;
}

bool MediumMagboltz::GetDeexcitationProduct(const unsigned int i, double& t, double& s,
                                            int& type, double& energy) const {

  if (i >= m_dxcProducts.size() || !(m_useDeexcitation || m_usePenning)) {
    return false;
  }
  t = m_dxcProducts[i].t;
  s = m_dxcProducts[i].s;
  type = m_dxcProducts[i].type;
  energy = m_dxcProducts[i].energy;
  return true;
}

bool MediumMagboltz::GetIonisationProduct(const unsigned int i, int& type,
                                          double& energy) const {

  if (i >= m_ionProducts.size()) {
    std::cerr << m_className << "::GetIonisationProduct:\n"
              << "    Index (" << i << ") out of range.\n";
    return false;
  }

  type = m_ionProducts[i].type;
  energy = m_ionProducts[i].energy;
  return true;
}

double MediumMagboltz::GetPhotonCollisionRate(const double e) {

  if (e <= 0.) {
    std::cerr << m_className << "::GetPhotonCollisionRate:\n";
    std::cerr << "    Photon energy must be greater than zero.\n";
    return m_cfTotGamma[0];
  }
  if (e > m_eFinalGamma && m_useAutoAdjust) {
    std::cerr << m_className << "::GetPhotonCollisionRate:\n";
    std::cerr << "    Collision rate at " << e
              << " eV is not included in the current table.\n";
    std::cerr << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxPhotonEnergy(1.05 * e);
  }

  if (m_isChanged) {
    if (!Mixer()) {
      std::cerr << m_className << "::GetPhotonCollisionRate:\n";
      std::cerr << "     Error calculating the collision rates table.\n";
      return 0.;
    }
    m_isChanged = false;
  }

  int iE = int(e / m_eStepGamma);
  if (iE >= nEnergyStepsGamma) iE = nEnergyStepsGamma - 1;
  if (iE < 0) iE = 0;

  double cfSum = m_cfTotGamma[iE];
  if (m_useDeexcitation && m_useRadTrap && !m_deexcitations.empty()) {
    // Loop over the excitations.
    const unsigned int nDeexcitations = m_deexcitations.size();
    for (unsigned int i = 0; i < nDeexcitations; ++i) {
      if (m_deexcitations[i].cf > 0. &&
          fabs(e - m_deexcitations[i].energy) <= m_deexcitations[i].width) {
        cfSum +=
            m_deexcitations[i].cf * TMath::Voigt(e - m_deexcitations[i].energy,
                                               m_deexcitations[i].sDoppler,
                                               2 * m_deexcitations[i].gPressure);
      }
    }
  }

  return cfSum;
}

bool MediumMagboltz::GetPhotonCollision(const double e, int& type, int& level,
                                        double& e1, double& ctheta, int& nsec,
                                        double& esec) {

  if (e > m_eFinalGamma && m_useAutoAdjust) {
    std::cerr << m_className << "::GetPhotonCollision:\n";
    std::cerr << "    Provided electron energy  (" << e
              << " eV) exceeds current energy range  (" << m_eFinalGamma
              << " eV).\n";
    std::cerr << "    Increasing energy range to " << 1.05 * e << " eV.\n";
    SetMaxPhotonEnergy(1.05 * e);
  } else if (e <= 0.) {
    std::cerr << m_className << "::GetPhotonCollision:\n";
    std::cerr << "    Photon energy must be greater than zero.\n";
    return false;
  }

  if (m_isChanged) {
    if (!Mixer()) {
      std::cerr << m_className << "::GetPhotonCollision:\n";
      std::cerr << "    Error calculating the collision rates table.\n";
      return false;
    }
    m_isChanged = false;
  }

  // Energy interval
  int iE = int(e / m_eStepGamma);
  if (iE >= nEnergyStepsGamma) iE = nEnergyStepsGamma - 1;
  if (iE < 0) iE = 0;

  double r = m_cfTotGamma[iE];
  if (m_useDeexcitation && m_useRadTrap && !m_deexcitations.empty()) {
    int nLines = 0;
    std::vector<double> pLine(0);
    std::vector<int> iLine(0);
    // Loop over the excitations.
    const unsigned int nDeexcitations = m_deexcitations.size();
    for (unsigned int i = 0; i < nDeexcitations; ++i) {
      if (m_deexcitations[i].cf > 0. &&
          fabs(e - m_deexcitations[i].energy) <= m_deexcitations[i].width) {
        r += m_deexcitations[i].cf * TMath::Voigt(e - m_deexcitations[i].energy,
                                                m_deexcitations[i].sDoppler,
                                                2 * m_deexcitations[i].gPressure);
        pLine.push_back(r);
        iLine.push_back(i);
        ++nLines;
      }
    }
    r *= RndmUniform();
    if (nLines > 0 && r >= m_cfTotGamma[iE]) {
      // Photon is absorbed by a discrete line.
      for (int i = 0; i < nLines; ++i) {
        if (r <= pLine[i]) {
          ++m_nPhotonCollisions[PhotonCollisionTypeExcitation];
          int fLevel = 0;
          ComputeDeexcitationInternal(iLine[i], fLevel);
          type = PhotonCollisionTypeExcitation;
          nsec = nDeexcitationProducts;
          return true;
        }
      }
      std::cerr << m_className << "::GetPhotonCollision:\n";
      std::cerr << "    Random sampling of deexcitation line failed.\n";
      std::cerr << "    Program bug!\n";
      return false;
    }
  } else {
    r *= RndmUniform();
  }

  int iLow = 0;
  int iUp = nPhotonTerms - 1;
  if (r <= m_cfGamma[iE][iLow]) {
    level = iLow;
  } else if (r >= m_cfGamma[iE][iUp]) {
    level = iUp;
  } else {
    int iMid;
    while (iUp - iLow > 1) {
      iMid = (iLow + iUp) >> 1;
      if (r < m_cfGamma[iE][iMid]) {
        iUp = iMid;
      } else {
        iLow = iMid;
      }
    }
    level = iUp;
  }

  nsec = 0;
  esec = e1 = 0.;
  type = csTypeGamma[level];
  // Collision type
  type = type % nCsTypesGamma;
  int ngas = int(csTypeGamma[level] / nCsTypesGamma);
  ++m_nPhotonCollisions[type];
  // Ionising collision
  if (type == 1) {
    esec = e - m_ionPot[ngas];
    if (esec < Small) esec = Small;
    nsec = 1;
  }

  // Determine the scattering angle
  ctheta = 2 * RndmUniform() - 1.;

  return true;
}

void MediumMagboltz::ResetCollisionCounters() {

  for (int j = nCsTypes; j--;) m_nCollisions[j] = 0;
  m_nCollisionsDetailed.resize(m_nTerms);
  for (unsigned int j = 0; j < m_nTerms; ++j) m_nCollisionsDetailed[j] = 0;
  m_nPenning = 0;
  for (int j = nCsTypesGamma; j--;) m_nPhotonCollisions[j] = 0;
}

unsigned int MediumMagboltz::GetNumberOfElectronCollisions() const {

  unsigned int ncoll = 0;
  for (int j = nCsTypes; j--;) ncoll += m_nCollisions[j];
  return ncoll;
}

unsigned int MediumMagboltz::GetNumberOfElectronCollisions(
    int& nElastic, int& nIonisation, int& nAttachment, int& nInelastic,
    int& nExcitation, int& nSuperelastic) const {

  nElastic = m_nCollisions[ElectronCollisionTypeElastic];
  nIonisation = m_nCollisions[ElectronCollisionTypeIonisation];
  nAttachment = m_nCollisions[ElectronCollisionTypeAttachment];
  nInelastic = m_nCollisions[ElectronCollisionTypeInelastic];
  nExcitation = m_nCollisions[ElectronCollisionTypeExcitation];
  nSuperelastic = m_nCollisions[ElectronCollisionTypeSuperelastic];
  return nElastic + nIonisation + nAttachment + nInelastic + nExcitation +
         nSuperelastic;
}

int MediumMagboltz::GetNumberOfLevels() {

  if (m_isChanged) {
    if (!Mixer()) {
      std::cerr << m_className << "::GetNumberOfLevels:\n";
      std::cerr << "    Error calculating the collision rates table.\n";
      return 0;
    }
    m_isChanged = false;
  }

  return m_nTerms;
}

bool MediumMagboltz::GetLevel(const unsigned int i, int& ngas, int& type,
                              std::string& descr, double& e) {

  if (m_isChanged) {
    if (!Mixer()) {
      std::cerr << m_className << "::GetLevel:\n";
      std::cerr << "    Error calculating the collision rates table.\n";
      return false;
    }
    m_isChanged = false;
  }

  if (i >= m_nTerms) {
    std::cerr << m_className << "::GetLevel:\n";
    std::cerr << "    Requested level (" << i << ") does not exist.\n";
    return false;
  }

  // Collision type
  type = m_csType[i] % nCsTypes;
  ngas = int(m_csType[i] / nCsTypes);
  // Description (from Magboltz)
  descr = std::string(50, ' ');
  for (int j = 50; j--;) descr[j] = m_description[i][j];
  // Threshold energy
  e = m_rgas[ngas] * m_energyLoss[i];
  if (m_debug) {
    std::cout << m_className << "::GetLevel:\n";
    std::cout << "    Level " << i << ": " << descr << "\n";
    std::cout << "    Type " << type << "\n",
        std::cout << "    Threshold energy: " << e << " eV\n";
    if (type == ElectronCollisionTypeExcitation && m_usePenning &&
        e > m_minIonPot) {
      std::cout << "    Penning transfer coefficient: " << m_rPenning[i] << "\n";
    } else if (type == ElectronCollisionTypeExcitation && m_useDeexcitation) {
      const int idxc = m_iDeexcitation[i];
      if (idxc < 0 || idxc >= (int)m_deexcitations.size()) {
        std::cout << "    Deexcitation cascade not implemented.\n";
        return true;
      }
      if (m_deexcitations[idxc].osc > 0.) {
        std::cout << "    Oscillator strength: " << m_deexcitations[idxc].osc
                  << "\n";
      }
      std::cout << "    Decay channels:\n";
      for (int j = 0; j < m_deexcitations[idxc].nChannels; ++j) {
        if (m_deexcitations[idxc].type[j] == DxcTypeRad) {
          std::cout << "      Radiative decay to ";
          if (m_deexcitations[idxc].final[j] < 0) {
            std::cout << "ground state: ";
          } else {
            std::cout << m_deexcitations[m_deexcitations[idxc].final[j]].label
                      << ": ";
          }
        } else if (m_deexcitations[idxc].type[j] == DxcTypeCollIon) {
          if (m_deexcitations[idxc].final[j] < 0) {
            std::cout << "      Penning ionisation: ";
          } else {
            std::cout << "      Associative ionisation: ";
          }
        } else if (m_deexcitations[idxc].type[j] == DxcTypeCollNonIon) {
          if (m_deexcitations[idxc].final[j] >= 0) {
            std::cout << "      Collision-induced transition to "
                      << m_deexcitations[m_deexcitations[idxc].final[j]].label
                      << ": ";
          } else {
            std::cout << "      Loss: ";
          }
        }
        if (j == 0) {
          std::cout << std::setprecision(5) << m_deexcitations[idxc].p[j] * 100.
                    << "%\n";
        } else {
          std::cout << std::setprecision(5) << (m_deexcitations[idxc].p[j] -
                                                m_deexcitations[idxc].p[j - 1]) *
                                                   100. << "%\n";
        }
      }
    }
  }

  return true;
}

unsigned int MediumMagboltz::GetNumberOfElectronCollisions(const unsigned int level) const {

  if (level >= m_nTerms) {
    std::cerr << m_className << "::GetNumberOfElectronCollisions:\n"
              << "    Cross-section term (" << level << ") does not exist.\n";
    return 0;
  }
  return m_nCollisionsDetailed[level];
}

int MediumMagboltz::GetNumberOfPhotonCollisions() const {

  int ncoll = 0;
  for (int j = nCsTypesGamma; j--;) ncoll += m_nPhotonCollisions[j];
  return ncoll;
}

int MediumMagboltz::GetNumberOfPhotonCollisions(int& nElastic, int& nIonising,
                                                int& nInelastic) const {

  nElastic = m_nPhotonCollisions[0];
  nIonising = m_nPhotonCollisions[1];
  nInelastic = m_nPhotonCollisions[2];
  return nElastic + nIonising + nInelastic;
}

bool MediumMagboltz::GetGasNumberMagboltz(const std::string& input,
                                          int& number) const {

  if (input == "") {
    number = 0;
    return false;
  }

  // CF4
  if (input == "CF4") {
    number = 1;
    return true;
  }
  // Argon
  if (input == "Ar") {
    number = 2;
    return true;
  }
  // Helium 4
  if (input == "He" || input == "He-4") {
    number = 3;
    return true;
  }
  // Helium 3
  if (input == "He-3") {
    number = 4;
    return true;
  }
  // Neon
  if (input == "Ne") {
    number = 5;
    return true;
  }
  // Krypton
  if (input == "Kr") {
    number = 6;
    return true;
  }
  // Xenon
  if (input == "Xe") {
    number = 7;
    return true;
  }
  // Methane
  if (input == "CH4") {
    number = 8;
    return true;
  }
  // Ethane
  if (input == "C2H6") {
    number = 9;
    return true;
  }
  // Propane
  if (input == "C3H8") {
    number = 10;
    return true;
  }
  // Isobutane
  if (input == "iC4H10") {
    number = 11;
    return true;
  }
  // Carbon dioxide (CO2)
  if (input == "CO2") {
    number = 12;
    return true;
  }
  // Neopentane
  if (input == "neoC5H12") {
    number = 13;
    return true;
  }
  // Water
  if (input == "H2O") {
    number = 14;
    return true;
  }
  // Oxygen
  if (input == "O2") {
    number = 15;
    return true;
  }
  // Nitrogen
  if (input == "N2") {
    number = 16;
    return true;
  }
  // Nitric oxide (NO)
  if (input == "NO") {
    number = 17;
    return true;
  }
  // Nitrous oxide (N2O)
  if (input == "N2O") {
    number = 18;
    return true;
  }
  // Ethene (C2H4)
  if (input == "C2H4") {
    number = 19;
    return true;
  }
  // Acetylene (C2H2)
  if (input == "C2H2") {
    number = 20;
    return true;
  }
  // Hydrogen
  if (input == "H2") {
    number = 21;
    return true;
  }
  // Deuterium
  if (input == "D2") {
    number = 22;
    return true;
  }
  // Carbon monoxide (CO)
  if (input == "CO") {
    number = 23;
    return true;
  }
  // Methylal (dimethoxymethane, CH3-O-CH2-O-CH3, "hot" version)
  if (input == "Methylal") {
    number = 24;
    return true;
  }
  // DME
  if (input == "DME") {
    number = 25;
    return true;
  }
  // Reid step
  if (input == "Reid-Step") {
    number = 26;
    return true;
  }
  // Maxwell model
  if (input == "Maxwell-Model") {
    number = 27;
    return true;
  }
  // Reid ramp
  if (input == "Reid-Ramp") {
    number = 28;
    return true;
  }
  // C2F6
  if (input == "C2F6") {
    number = 29;
    return true;
  }
  // SF6
  if (input == "SF6") {
    number = 30;
    return true;
  }
  // NH3
  if (input == "NH3") {
    number = 31;
    return true;
  }
  // Propene
  if (input == "C3H6") {
    number = 32;
    return true;
  }
  // Cyclopropane
  if (input == "cC3H6") {
    number = 33;
    return true;
  }
  // Methanol
  if (input == "CH3OH") {
    number = 34;
    return true;
  }
  // Ethanol
  if (input == "C2H5OH") {
    number = 35;
    return true;
  }
  // Propanol
  if (input == "C3H7OH") {
    number = 36;
    return true;
  }
  // Cesium / Caesium.
  if (input == "Cs") {
    number = 37;
    return true;
  }
  // Fluorine
  if (input == "F2") {
    number = 38;
    return true;
  }
  if (input == "CS2") {
    number = 39;
    return true;
  }
  // COS
  if (input == "COS") {
    number = 40;
    return true;
  }
  // Deuterated methane
  if (input == "CD4") {
    number = 41;
    return true;
  }
  // BF3
  if (input == "BF3") {
    number = 42;
    return true;
  }
  // C2H2F4 (C2HF5).
  if (input == "C2HF5" || input == "C2H2F4") {
    number = 43;
    return true;
  }
  // TMA
  if (input == "TMA") {
    number = 44;
    return true;
  }
  // CHF3
  if (input == "CHF3") {
    number = 50;
    return true;
  }
  // CF3Br
  if (input == "CF3Br") {
    number = 51;
    return true;
  }
  // C3F8
  if (input == "C3F8") {
    number = 52;
    return true;
  }
  // Ozone
  if (input == "O3") {
    number = 53;
    return true;
  }
  // Mercury
  if (input == "Hg") {
    number = 54;
    return true;
  }
  // H2S
  if (input == "H2S") {
    number = 55;
    return true;
  }
  // n-Butane
  if (input == "nC4H10") {
    number = 56;
    return true;
  }
  // n-Pentane
  if (input == "nC5H12") {
    number = 57;
    return true;
  }
  // Nitrogen
  if (input == "N2 (Phelps)") {
    number = 58;
    return true;
  }
  // Germane, GeH4
  if (input == "GeH4") {
    number = 59;
    return true;
  }
  // Silane, SiH4
  if (input == "SiH4") {
    number = 60;
    return true;
  }

  std::cerr << m_className << "::GetGasNumberMagboltz:\n";
  std::cerr << "    Gas " << input << " is not defined.\n";
  return false;
}

bool MediumMagboltz::Mixer(const bool verbose) {

  // Set constants and parameters in Magboltz common blocks.
  Magboltz::cnsts_.echarg = ElementaryCharge * 1.e-15;
  Magboltz::cnsts_.emass = ElectronMassGramme;
  Magboltz::cnsts_.amu = AtomicMassUnit;
  Magboltz::cnsts_.pir2 = BohrRadius * BohrRadius * Pi;
  Magboltz::inpt_.ary = RydbergEnergy;

  Magboltz::inpt_.akt = BoltzmannConstant * m_temperature;
  Magboltz::inpt_.tempc = m_temperature - ZeroCelsius;
  Magboltz::inpt_.torr = m_pressure;

  Magboltz::inpt_.nGas = m_nComponents;
  Magboltz::inpt_.nStep = nEnergySteps;
  if (m_useAnisotropic) {
    Magboltz::inpt_.nAniso = 2;
  } else {
    Magboltz::inpt_.nAniso = 0;
  }

  // Calculate the atomic density (ideal gas law).
  const double dens = GetNumberDensity();
  // Prefactor for calculation of scattering rate from cross-section.
  const double prefactor = dens * SpeedOfLight * sqrt(2. / ElectronMass);

  // Fill the electron energy array, reset the collision rates.
  for (int i = nEnergySteps; i--;) {
    m_cfTot[i] = 0.;
    for (int j = nMaxLevels; j--;) {
      m_cf[i][j] = 0.;
      m_scatParameter[i][j] = 0.5;
      m_scatCut[i][j] = 1.;
    }
  }
  for (int i = nEnergyStepsLog; i--;) {
    m_cfTotLog[i] = 0.;
    for (int j = nMaxLevels; j--;) {
      m_cfLog[i][j] = 0.;
      m_scatParameter[i][j] = 0.5;
      m_scatCut[i][j] = 1.;
    }
  }

  m_deexcitations.clear();
  for (int i = nMaxLevels; i--;) {
    m_scatModel[i] = 0;
    m_iDeexcitation[i] = -1;
    m_wOpalBeaty[i] = 1.;
  }

  m_minIonPot = -1.;
  for (unsigned int i = 0; i < m_nMaxGases; ++i) {
    m_ionPot[i] = -1.;
    m_gsGreenSawada[i] = 1.;
    m_gbGreenSawada[i] = 0.;
    m_tsGreenSawada[i] = 0.;
    m_taGreenSawada[i] = 0.;
    m_tbGreenSawada[i] = 0.;
    m_hasGreenSawada[i] = false;
  }
  // Cross-sections
  // 0: total, 1: elastic,
  // 2: ionisation, 3: attachment,
  // 4, 5: unused
  static double q[nEnergySteps][6];
  // Parameters for scattering angular distribution
  static double pEqEl[nEnergySteps][6];
  // Inelastic cross-sections
  static double qIn[nEnergySteps][nMaxInelasticTerms];
  // Ionisation cross-sections
  static double qIon[nEnergySteps][8];
  // Parameters for angular distribution in inelastic collisions
  static double pEqIn[nEnergySteps][nMaxInelasticTerms];
  // Parameters for angular distribution in ionising collisions
  static double pEqIon[nEnergySteps][8];
  // Opal-Beaty parameter
  static double eoby[nEnergySteps];
  // Penning transfer parameters
  static double penFra[nMaxInelasticTerms][3];
  // Description of cross-section terms
  static char scrpt[260][50];

  // Check the gas composition and establish the gas numbers.
  int gasNumber[m_nMaxGases];
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (!GetGasNumberMagboltz(m_gas[i], gasNumber[i])) {
      std::cerr << m_className << "::Mixer:\n";
      std::cerr << "    Gas " << m_gas[i] << " has no corresponding"
                << " gas number in Magboltz.\n";
      return false;
    }
  }

  if (m_debug || verbose) {
    std::cout << m_className << "::Mixer:\n";
    std::cout << "    Creating table of collision rates with\n";
    std::cout << "    " << nEnergySteps << " linear energy steps between 0 and "
              << std::min(m_eFinal, m_eHigh) << " eV\n";
    if (m_eFinal > m_eHigh) {
      std::cout << "    " << nEnergyStepsLog
                << " logarithmic energy steps between " << m_eHigh << " and "
                << m_eFinal << " eV\n";
    }
  }
  m_nTerms = 0;

  std::ofstream outfile;
  if (m_useCsOutput) {
    outfile.open("cs.txt", std::ios::out);
    outfile << "# energy [eV] vs. cross-section [cm2]\n";
  }

  // Loop over the gases in the mixture.
  for (unsigned int iGas = 0; iGas < m_nComponents; ++iGas) {
    if (m_eFinal <= m_eHigh) {
      Magboltz::inpt_.efinal = m_eFinal;
    } else {
      Magboltz::inpt_.efinal = m_eHigh;
    }
    Magboltz::inpt_.estep = m_eStep;

    // Number of inelastic cross-section terms
    long long nIn = 0;
    long long nIon = 0;
    // Threshold energies
    double e[6] = {0., 0., 0., 0., 0., 0.};
    double eIn[nMaxInelasticTerms] = {0.};
    double eIon[8] = {0.};
    // Virial coefficient (not used)
    double virial = 0.;
    // Scattering algorithms
    long long kIn[nMaxInelasticTerms] = {0};
    long long kEl[6] = {0, 0, 0, 0, 0, 0};
    char name[] = "                         ";

    // Retrieve the cross-section data for this gas from Magboltz.
    long long ngs = gasNumber[iGas];
    Magboltz::gasmix_(&ngs, q[0], qIn[0], &nIn, e, eIn, name, &virial, eoby,
                      pEqEl[0], pEqIn[0], penFra[0], kEl, kIn, qIon[0],
                      pEqIon[0], eIon, &nIon, scrpt);
    if (m_debug || verbose) {
      const double massAmu =
          (2. / e[1]) * ElectronMass / AtomicMassUnitElectronVolt;
      std::cout << "    " << name << "\n";
      std::cout << "      mass:                 " << massAmu << " amu\n";
      if (nIon > 1) {
        std::cout << "      ionisation threshold: " << eIon[0] << " eV\n";
      } else {
        std::cout << "      ionisation threshold: " << e[2] << " eV\n";
      }
      if (e[3] > 0. && e[4] > 0.) {
        std::cout << "      cross-sections at minimum ionising energy:\n";
        std::cout << "        excitation: " << e[3] * 1.e18 << " Mbarn\n";
        std::cout << "        ionisation: " << e[4] * 1.e18 << " Mbarn\n";
      }
    }
    int np0 = m_nTerms;

    // Make sure there is still sufficient space.
    if (np0 + nIn + nIon + 1 >= nMaxLevels) {
      std::cerr << m_className << "::Mixer:\n";
      std::cerr << "    Max. number of levels (" << nMaxLevels
                << ") exceeded.\n";
      return false;
    }

    double van = m_fraction[iGas] * prefactor;

    int np = np0;
    if (m_useCsOutput) {
      outfile << "# cross-sections for " << name << "\n";
      outfile << "# cross-section types:\n";
      outfile << "# elastic\n";
    }
    // Elastic scattering
    ++m_nTerms;
    m_scatModel[np] = kEl[1];
    const double r = 1. + e[1] / 2.;
    m_rgas[iGas] = r;
    m_energyLoss[np] = 0.;
    for (int j = 0; j < 50; ++j) {
      m_description[np][j] = scrpt[1][j];
    }
    m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeElastic;
    bool withIon = false;
    // Ionisation
    if (nIon > 1) {
      for (int j = 0; j < nIon; ++j) {
        if (m_eFinal < eIon[j]) continue;
        withIon = true;
        ++m_nTerms;
        ++np;
        m_scatModel[np] = kEl[2];
        m_energyLoss[np] = eIon[j] / r;
        // TODO
        m_wOpalBeaty[np] = eoby[j];
        if (m_gas[iGas] == "CH4") {
          if (fabs(eIon[j] - 21.) < 0.1) {
            m_wOpalBeaty[np] = 14.;
          } else if (fabs(eIon[j] - 291.) < 0.1) {
            m_wOpalBeaty[np] = 200.;
          }
        }
        for (int k = 0; k < 50; ++k) {
          m_description[np][k] = scrpt[2 + j][k];
        }
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeIonisation;
        if (m_useCsOutput) {
          outfile << "# " << m_description[np] << "\n";
        }
      }
      m_gsGreenSawada[iGas] = eoby[0];
      m_tbGreenSawada[iGas] = 2 * eIon[0];
      m_ionPot[iGas] = eIon[0];
    } else {
      if (m_eFinal >= e[2]) {
        withIon = true;
        ++m_nTerms;
        ++np;
        m_scatModel[np] = kEl[2];
        m_energyLoss[np] = e[2] / r;
        m_wOpalBeaty[np] = eoby[0];
        m_gsGreenSawada[iGas] = eoby[0];
        m_tbGreenSawada[iGas] = 2 * e[2];
        m_ionPot[iGas] = e[2];
        for (int j = 0; j < 50; ++j) {
          m_description[np][j] = scrpt[2][j];
        }
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeIonisation;
        if (m_useCsOutput) {
          outfile << "# ionisation (gross)\n";
        }
      }
    }
    // Attachment
    ++m_nTerms;
    ++np;
    m_scatModel[np] = 0;
    m_energyLoss[np] = 0.;
    for (int j = 0; j < 50; ++j) {
      m_description[np][j] = scrpt[2 + nIon][j];
    }
    m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeAttachment;
    if (m_useCsOutput) {
      outfile << "# attachment\n";
    }
    // Inelastic terms
    int nExc = 0, nSuperEl = 0;
    for (int j = 0; j < nIn; ++j) {
      ++np;
      m_scatModel[np] = kIn[j];
      m_energyLoss[np] = eIn[j] / r;
      for (int k = 0; k < 50; ++k) {
        m_description[np][k] = scrpt[5 + nIon + j][k];
      }
      if ((m_description[np][1] == 'E' && m_description[np][2] == 'X') ||
          (m_description[np][0] == 'E' && m_description[np][1] == 'X') ||
          (m_gas[iGas] == "N2" && eIn[j] > 6.)) {
        // Excitation
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeExcitation;
        ++nExc;
      } else if (eIn[j] < 0.) {
        // Super-elastic collision
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeSuperelastic;
        ++nSuperEl;
      } else {
        // Inelastic collision
        m_csType[np] = nCsTypes * iGas + ElectronCollisionTypeInelastic;
      }
      if (m_useCsOutput) {
        outfile << "# " << m_description[np] << "\n";
      }
    }
    m_nTerms += nIn;
    // Loop over the energy table.
    for (int iE = 0; iE < nEnergySteps; ++iE) {
      np = np0;
      if (m_useCsOutput) {
        outfile << (iE + 0.5) * m_eStep << "  " << q[iE][1] << "  ";
      }
      // Elastic scattering
      m_cf[iE][np] = q[iE][1] * van;
      if (m_scatModel[np] == 1) {
        ComputeAngularCut(pEqEl[iE][1], m_scatCut[iE][np], m_scatParameter[iE][np]);
      } else if (m_scatModel[np] == 2) {
        m_scatParameter[iE][np] = pEqEl[iE][1];
      }
      // Ionisation
      if (withIon) {
        if (nIon > 1) {
          for (int j = 0; j < nIon; ++j) {
            if (m_eFinal < eIon[j]) continue;
            ++np;
            m_cf[iE][np] = qIon[iE][j] * van;
            if (m_scatModel[np] == 1) {
              ComputeAngularCut(pEqIon[iE][j], m_scatCut[iE][np],
                                m_scatParameter[iE][np]);
            } else if (m_scatModel[np] == 2) {
              m_scatParameter[iE][np] = pEqIon[iE][j];
            }
            if (m_useCsOutput) {
              outfile << qIon[iE][j] << "  ";
            }
          }
        } else {
          ++np;
          m_cf[iE][np] = q[iE][2] * van;
          if (m_scatModel[np] == 1) {
            ComputeAngularCut(pEqEl[iE][2], m_scatCut[iE][np],
                              m_scatParameter[iE][np]);
          } else if (m_scatModel[np] == 2) {
            m_scatParameter[iE][np] = pEqEl[iE][2];
          }
          if (m_useCsOutput) {
            outfile << q[iE][2] << "  ";
          }
        }
      }
      // Attachment
      ++np;
      m_cf[iE][np] = q[iE][3] * van;
      m_scatParameter[iE][np] = 0.5;
      if (m_useCsOutput) {
        outfile << q[iE][3] << "  ";
      }
      // Inelastic terms
      for (int j = 0; j < nIn; ++j) {
        ++np;
        if (m_useCsOutput) outfile << qIn[iE][j] << "  ";
        m_cf[iE][np] = qIn[iE][j] * van;
        // Scale the excitation cross-sections (for error estimates).
        m_cf[iE][np] *= m_scaleExc[iGas];
        // Temporary hack for methane dissociative excitations:
        if (m_description[np][5] == 'D' && m_description[np][6] == 'I' &&
            m_description[np][7] == 'S') {
          // if ((iE + 0.5) * m_eStep > 40.) {
          //   m_cf[iE][np] *= 0.8;
          // } else if ((iE + 0.5) * m_eStep > 30.) {
          //   m_cf[iE][np] *= (1. - ((iE + 0.5) * m_eStep - 30.) * 0.02);
          // }
        }
        if (m_cf[iE][np] < 0.) {
          std::cerr << m_className << "::Mixer:\n";
          std::cerr << "    Negative inelastic cross-section at "
                    << (iE + 0.5) * m_eStep << " eV.\n";
          std::cerr << "    Set to zero.\n";
          m_cf[iE][np] = 0.;
        }
        if (m_scatModel[np] == 1) {
          ComputeAngularCut(pEqIn[iE][j], m_scatCut[iE][np],
                            m_scatParameter[iE][np]);
        } else if (m_scatModel[np] == 2) {
          m_scatParameter[iE][np] = pEqIn[iE][j];
        }
      }
      if ((m_debug || verbose) && nIn > 0 && iE == nEnergySteps - 1) {
        std::cout << "      " << nIn << " inelastic terms (" << nExc
                  << " excitations, " << nSuperEl << " superelastic, "
                  << nIn - nExc - nSuperEl << " other)\n";
      }
      if (m_useCsOutput) outfile << "\n";
    }

    if (m_eFinal <= m_eHigh) continue;
    // Fill the high-energy part (logarithmic binning).
    // Calculate the growth factor.
    const double rLog = pow(m_eFinal / m_eHigh, 1. / nEnergyStepsLog);
    m_lnStep = log(rLog);
    // Set the upper limit of the first bin.
    double emax = m_eHigh * rLog;
    int imax = nEnergySteps - 1;
    for (int iE = 0; iE < nEnergyStepsLog; ++iE) {
      Magboltz::inpt_.estep = emax / (nEnergySteps - 0.5);
      Magboltz::inpt_.efinal = emax + 0.5 * Magboltz::inpt_.estep;
      Magboltz::gasmix_(&ngs, q[0], qIn[0], &nIn, e, eIn, name, &virial, eoby,
                        pEqEl[0], pEqIn[0], penFra[0], kEl, kIn, qIon[0],
                        pEqIon[0], eIon, &nIon, scrpt);
      np = np0;
      if (m_useCsOutput) {
        outfile << emax << "  " << q[imax][1] << "  ";
      }
      // Elastic scattering
      m_cfLog[iE][np] = q[imax][1] * van;
      if (m_scatModel[np] == 1) {
        ComputeAngularCut(pEqEl[imax][1], m_scatCutLog[iE][np],
                          m_scatParameterLog[iE][np]);
      } else if (m_scatModel[np] == 2) {
        m_scatParameterLog[iE][np] = pEqEl[imax][1];
      }
      // Ionisation
      if (withIon) {
        if (nIon > 1) {
          for (int j = 0; j < nIon; ++j) {
            if (m_eFinal < eIon[j]) continue;
            ++np;
            m_cfLog[iE][np] = qIon[imax][j] * van;
            if (m_scatModel[np] == 1) {
              ComputeAngularCut(pEqIon[imax][j], m_scatCutLog[iE][np],
                                m_scatParameterLog[iE][np]);
            } else if (m_scatModel[np] == 2) {
              m_scatParameterLog[iE][np] = pEqIon[imax][j];
            }
            if (m_useCsOutput) {
              outfile << qIon[imax][j] << "  ";
            }
          }
        } else {
          ++np;
          // Gross cross-section
          m_cfLog[iE][np] = q[imax][2] * van;
          // Counting cross-section
          // m_cfLog[iE][np] = q[imax][4] * van;
          if (m_scatModel[np] == 1) {
            ComputeAngularCut(pEqEl[imax][2], m_scatCutLog[iE][np],
                              m_scatParameterLog[iE][np]);
          } else if (m_scatModel[np] == 2) {
            m_scatParameterLog[iE][np] = pEqEl[imax][2];
          }
        }
      }
      // Attachment
      ++np;
      m_cfLog[iE][np] = q[imax][3] * van;
      if (m_useCsOutput) {
        outfile << q[imax][3] << "  ";
      }
      // Inelastic terms
      for (int j = 0; j < nIn; ++j) {
        ++np;
        if (m_useCsOutput) outfile << qIn[imax][j] << "  ";
        m_cfLog[iE][np] = qIn[imax][j] * van;
        // Scale the excitation cross-sections (for error estimates).
        m_cfLog[iE][np] *= m_scaleExc[iGas];
        if (m_cfLog[iE][np] < 0.) {
          std::cerr << m_className << "::Mixer:\n";
          std::cerr << "    Negative inelastic cross-section at " << emax
                    << " eV.\n";
          std::cerr << "    Set to zero.\n";
          m_cfLog[iE][np] = 0.;
        }
        if (m_scatModel[np] == 1) {
          ComputeAngularCut(pEqIn[imax][j], m_scatCutLog[iE][np],
                            m_scatParameterLog[iE][np]);
        } else if (m_scatModel[np] == 2) {
          m_scatParameterLog[iE][np] = pEqIn[imax][j];
        }
      }
      if (m_useCsOutput) outfile << "\n";
      // Increase the energy.
      emax *= rLog;
    }
  }
  if (m_useCsOutput) outfile.close();

  // Find the smallest ionisation threshold.
  std::string minIonPotGas = "";
  for (unsigned int i = 0; i < m_nMaxGases; ++i) {
    if (m_ionPot[i] < 0.) continue;
    if (m_minIonPot < 0.) {
      m_minIonPot = m_ionPot[i];
      minIonPotGas = m_gas[i];
    } else if (m_ionPot[i] < m_minIonPot) {
      m_minIonPot = m_ionPot[i];
      minIonPotGas = m_gas[i];
    }
  }

  if (m_debug || verbose) {
    std::cout << m_className << "::Mixer:\n";
    std::cout << "    Lowest ionisation threshold in the mixture:\n";
    std::cout << "      " << m_minIonPot << " eV (" << minIonPotGas << ")\n";
  }

  for (int iE = nEnergySteps; iE--;) {
    // Calculate the total collision frequency.
    for (unsigned int k = 0; k < m_nTerms; ++k) {
      if (m_cf[iE][k] < 0.) {
        std::cerr << m_className << "::Mixer:\n";
        std::cerr << "    Negative collision rate at " << (iE + 0.5) * m_eStep
                  << " eV. Set to zero.\n";
        m_cf[iE][k] = 0.;
      }
      m_cfTot[iE] += m_cf[iE][k];
    }
    // Normalise the collision probabilities.
    if (m_cfTot[iE] > 0.) {
      for (unsigned int k = 0; k < m_nTerms; ++k) m_cf[iE][k] /= m_cfTot[iE];
    }
    for (unsigned int k = 1; k < m_nTerms; ++k) {
      m_cf[iE][k] += m_cf[iE][k - 1];
    }
    const double ekin = m_eStep * (iE + 0.5);
    m_cfTot[iE] *= sqrt(ekin);
    // Use relativistic expression at high energies.
    if (ekin > 1.e3) {
      m_cfTot[iE] *=
          sqrt(1. + 0.5 * ekin / ElectronMass) / (1. + ekin / ElectronMass);
    }
  }

  if (m_eFinal > m_eHigh) {
    const double rLog = pow(m_eFinal / m_eHigh, 1. / nEnergyStepsLog);
    for (int iE = nEnergyStepsLog; iE--;) {
      // Calculate the total collision frequency.
      for (unsigned int k = 0; k < m_nTerms; ++k) {
        if (m_cfLog[iE][k] < 0.) m_cfLog[iE][k] = 0.;
        m_cfTotLog[iE] += m_cfLog[iE][k];
      }
      // Normalise the collision probabilities.
      if (m_cfTotLog[iE] > 0.) {
        for (int k = m_nTerms; k--;) m_cfLog[iE][k] /= m_cfTotLog[iE];
      }
      for (unsigned int k = 1; k < m_nTerms; ++k) {
        m_cfLog[iE][k] += m_cfLog[iE][k - 1];
      }
      const double ekin = m_eHigh * pow(rLog, iE + 1);
      m_cfTotLog[iE] *= sqrt(ekin) * sqrt(1. + 0.5 * ekin / ElectronMass) /
                      (1. + ekin / ElectronMass);
      // Store the logarithm (for log-log interpolation)
      m_cfTotLog[iE] = log(m_cfTotLog[iE]);
    }
  }

  // Determine the null collision frequency.
  m_cfNull = 0.;
  for (int j = 0; j < nEnergySteps; ++j) {
    if (m_cfTot[j] > m_cfNull) m_cfNull = m_cfTot[j];
  }
  if (m_eFinal > m_eHigh) {
    for (int j = 0; j < nEnergyStepsLog; ++j) {
      const double r = exp(m_cfTotLog[j]);
      if (r > m_cfNull) m_cfNull = r;
    }
  }

  // Reset the collision counters.
  m_nCollisionsDetailed.resize(m_nTerms);
  for (int j = nCsTypes; j--;) m_nCollisions[j] = 0;
  for (int j = m_nTerms; j--;) m_nCollisionsDetailed[j] = 0;

  if (m_debug || verbose) {
    std::cout << m_className << "::Mixer:\n";
    std::cout << "    Energy [eV]    Collision Rate [ns-1]\n";
    for (int i = 0; i < 8; ++i) {
      const double emax = std::min(m_eHigh, m_eFinal);
      std::cout << "    " << std::fixed << std::setw(10) << std::setprecision(2)
                << (2 * i + 1) * emax / 16 << "    " << std::setw(18)
                << std::setprecision(2) << m_cfTot[(i + 1) * nEnergySteps / 16]
                << "\n";
    }
    std::cout << std::resetiosflags(std::ios_base::floatfield);
  }

  // Set up the de-excitation channels.
  if (m_useDeexcitation) {
    ComputeDeexcitationTable(verbose);
    const unsigned int nDeexcitations = m_deexcitations.size();
    for (unsigned int j = 0; j < nDeexcitations; ++j) {
      const int probCount = m_deexcitations[j].p.size();
      const int flvlCount = m_deexcitations[j].final.size();
      const int typeCount = m_deexcitations[j].type.size();
      if (!(probCount == flvlCount && flvlCount == typeCount &&
            typeCount == m_deexcitations[j].nChannels)) {
        std::cerr << m_className << "::Mixer:\n";
        std::cerr << "    Mismatch in deexcitation channel count.\n";
        std::cerr << "    Program bug!\n";
        std::cerr << "    Deexcitation handling is switched off.\n";
        m_useDeexcitation = false;
      }
    }
  }

  // Fill the photon collision rates table.
  if (!ComputePhotonCollisionTable(verbose)) {
    std::cerr << m_className << "::Mixer:\n";
    std::cerr << "    Photon collision rates could not be calculated.\n";
    if (m_useDeexcitation) {
      std::cerr << "    Deexcitation handling is switched off.\n";
      m_useDeexcitation = false;
    }
  }

  // Reset the Penning transfer parameters.
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    m_rPenning[i] = m_rPenningGlobal;
    int iGas = int(m_csType[i] / nCsTypes);
    if (m_rPenningGas[iGas] > Small) {
      m_rPenning[i] = m_rPenningGas[iGas];
      m_lambdaPenning[i] = m_lambdaPenningGas[iGas];
    }
  }

  // Set the Green-Sawada splitting function parameters.
  SetupGreenSawada();

  return true;
}

void MediumMagboltz::SetupGreenSawada() {

  for (unsigned int i = 0; i < m_nComponents; ++i) {
    m_taGreenSawada[i] = 1000.;
    m_hasGreenSawada[i] = true;
    if (m_gas[i] == "He" || m_gas[i] == "He-3") {
      m_tsGreenSawada[i] = -2.25;
      m_gsGreenSawada[i] = 15.5;
      m_gbGreenSawada[i] = 24.5;
    } else if (m_gas[i] == "Ne") {
      m_tsGreenSawada[i] = -6.49;
      m_gsGreenSawada[i] = 24.3;
      m_gbGreenSawada[i] = 21.6;
    } else if (m_gas[i] == "Ar") {
      m_tsGreenSawada[i] = 6.87;
      m_gsGreenSawada[i] = 6.92;
      m_gbGreenSawada[i] = 7.85;
    } else if (m_gas[i] == "Kr") {
      m_tsGreenSawada[i] = 3.90;
      m_gsGreenSawada[i] = 7.95;
      m_gbGreenSawada[i] = 13.5;
    } else if (m_gas[i] == "Xe") {
      m_tsGreenSawada[i] = 3.81;
      m_gsGreenSawada[i] = 7.93;
      m_gbGreenSawada[i] = 11.5;
    } else if (m_gas[i] == "H2" || m_gas[i] == "D2") {
      m_tsGreenSawada[i] = 1.87;
      m_gsGreenSawada[i] = 7.07;
      m_gbGreenSawada[i] = 7.7;
    } else if (m_gas[i] == "N2") {
      m_tsGreenSawada[i] = 4.71;
      m_gsGreenSawada[i] = 13.8;
      m_gbGreenSawada[i] = 15.6;
    } else if (m_gas[i] == "O2") {
      m_tsGreenSawada[i] = 1.86;
      m_gsGreenSawada[i] = 18.5;
      m_gbGreenSawada[i] = 12.1;
    } else if (m_gas[i] == "CH4") {
      m_tsGreenSawada[i] = 3.45;
      m_gsGreenSawada[i] = 7.06;
      m_gbGreenSawada[i] = 12.5;
    } else if (m_gas[i] == "H20") {
      m_tsGreenSawada[i] = 1.28;
      m_gsGreenSawada[i] = 12.8;
      m_gbGreenSawada[i] = 12.6;
    } else if (m_gas[i] == "CO") {
      m_tsGreenSawada[i] = 2.03;
      m_gsGreenSawada[i] = 13.3;
      m_gbGreenSawada[i] = 14.0;
    } else if (m_gas[i] == "C2H2") {
      m_tsGreenSawada[i] = 1.37;
      m_gsGreenSawada[i] = 9.28;
      m_gbGreenSawada[i] = 5.8;
    } else if (m_gas[i] == "NO") {
      m_tsGreenSawada[i] = -4.30;
      m_gsGreenSawada[i] = 10.4;
      m_gbGreenSawada[i] = 9.5;
    } else if (m_gas[i] == "CO2") {
      m_tsGreenSawada[i] = -2.46;
      m_gsGreenSawada[i] = 12.3;
      m_gbGreenSawada[i] = 13.8;
    } else {
      m_taGreenSawada[i] = 0.;
      m_hasGreenSawada[i] = false;
      if (m_useGreenSawada) {
        std::cout << m_className << "::SetupGreenSawada:\n"
                  << "    Fit parameters for " << m_gas[i] << " not available.\n"
                  << "    Opal-Beaty formula is used instead.\n";
      }
    }
  }
}

void MediumMagboltz::ComputeAngularCut(const double parIn, double& cut,
                                       double& parOut) const {

  // Set cuts on angular distribution and
  // renormalise forward scattering probability.

  if (parIn <= 1.) {
    cut = 1.;
    parOut = parIn;
    return;
  }

  const double rads = 2. / Pi;
  const double cns = parIn - 0.5;
  const double thetac = asin(2. * sqrt(cns - cns * cns));
  const double fac = (1. - cos(thetac)) / pow(sin(thetac), 2.);
  parOut = cns * fac + 0.5;
  cut = thetac * rads;
}

void MediumMagboltz::ComputeDeexcitationTable(const bool verbose) {

  for (int i = nMaxLevels; i--;) m_iDeexcitation[i] = -1;
  m_deexcitations.clear();

  // Optical data (for quencher photoabsorption cs and ionisation yield)
  OpticalData optData;

  // Presence flags, concentrations and indices of "de-excitable" gases.
  bool withAr = false;
  int iAr = 0;
  double cAr = 0.;
  bool withNe = false;

  std::map<std::string, int> mapLevels;
  // Make a mapping of all excitation levels.
  for (unsigned int i = 0; i < m_nTerms; ++i) {
    // Check if the level is an excitation.
    if (m_csType[i] % nCsTypes != ElectronCollisionTypeExcitation) continue;
    // Extract the index of the gas.
    const int ngas = int(m_csType[i] / nCsTypes);
    if (m_gas[ngas] == "Ar") {
      // Argon
      if (!withAr) {
        withAr = true;
        iAr = ngas;
        cAr = m_fraction[iAr];
      }
      // Get the level description (as specified in Magboltz).
      std::string level = "       ";
      for (int j = 0; j < 7; ++j) level[j] = m_description[i][5 + j];
      if (level == "1S5    ") mapLevels["Ar_1S5"] = i;
      else if (level == "1S4    ") mapLevels["Ar_1S4"] = i;
      else if (level == "1S3    ") mapLevels["Ar_1S3"] = i;
      else if (level == "1S2    ") mapLevels["Ar_1S2"] = i;
      else if (level == "2P10   ") mapLevels["Ar_2P10"] = i;
      else if (level == "2P9    ") mapLevels["Ar_2P9"] = i;
      else if (level == "2P8    ") mapLevels["Ar_2P8"] = i;
      else if (level == "2P7    ") mapLevels["Ar_2P7"] = i;
      else if (level == "2P6    ") mapLevels["Ar_2P6"] = i;
      else if (level == "2P5    ") mapLevels["Ar_2P5"] = i;
      else if (level == "2P4    ") mapLevels["Ar_2P4"] = i;
      else if (level == "2P3    ") mapLevels["Ar_2P3"] = i;
      else if (level == "2P2    ") mapLevels["Ar_2P2"] = i;
      else if (level == "2P1    ") mapLevels["Ar_2P1"] = i;
      else if (level == "3D6    ") mapLevels["Ar_3D6"] = i;
      else if (level == "3D5    ") mapLevels["Ar_3D5"] = i;
      else if (level == "3D3    ") mapLevels["Ar_3D3"] = i;
      else if (level == "3D4!   ") mapLevels["Ar_3D4!"] = i;
      else if (level == "3D4    ") mapLevels["Ar_3D4"] = i;
      else if (level == "3D1!!  ") mapLevels["Ar_3D1!!"] = i;
      else if (level == "2S5    ") mapLevels["Ar_2S5"] = i;
      else if (level == "2S4    ") mapLevels["Ar_2S4"] = i;
      else if (level == "3D1!   ") mapLevels["Ar_3D1!"] = i;
      else if (level == "3D2    ") mapLevels["Ar_3D2"] = i;
      else if (level == "3S1!!!!") mapLevels["Ar_3S1!!!!"] = i;
      else if (level == "3S1!!  ") mapLevels["Ar_3S1!!"] = i;
      else if (level == "3S1!!! ") mapLevels["Ar_3S1!!!"] = i;
      else if (level == "2S3    ") mapLevels["Ar_2S3"] = i;
      else if (level == "2S2    ") mapLevels["Ar_2S2"] = i;
      else if (level == "3S1!   ") mapLevels["Ar_3S1!"] = i;
      else if (level == "4D5    ") mapLevels["Ar_4D5"] = i;
      else if (level == "3S4    ") mapLevels["Ar_3S4"] = i;
      else if (level == "4D2    ") mapLevels["Ar_4D2"] = i;
      else if (level == "4S1!   ") mapLevels["Ar_4S1!"] = i;
      else if (level == "3S2    ") mapLevels["Ar_3S2"] = i;
      else if (level == "5D5    ") mapLevels["Ar_5D5"] = i;
      else if (level == "4S4    ") mapLevels["Ar_4S4"] = i;
      else if (level == "5D2    ") mapLevels["Ar_5D2"] = i;
      else if (level == "6D5    ") mapLevels["Ar_6D5"] = i;
      else if (level == "5S1!   ") mapLevels["Ar_5S1!"] = i;
      else if (level == "4S2    ") mapLevels["Ar_4S2"] = i;
      else if (level == "5S4    ") mapLevels["Ar_5S4"] = i;
      else if (level == "6D2    ") mapLevels["Ar_6D2"] = i;
      else if (level == "HIGH   ") mapLevels["Ar_Higher"] = i;
      else {
        std::cerr << m_className << "::ComputeDeexcitationTable:\n";
        std::cerr << "    Unknown excitation level:\n";
        std::cerr << "      Ar " << level << "\n";
      }
    } else if (m_gas[ngas] == "Ne") {
      // Neon
      if (!withNe) {
        withNe = true;
      }
      std::string level = "       ";
      for (int j = 0; j < 7; ++j) level[j] = m_description[i][3 + j];
      if (level == "  1S5  ") mapLevels["Ne_1S5"] = i;
      else if (level == "  1S4  ") mapLevels["Ne_1S4"] = i;
      else if (level == "  1S3  ") mapLevels["Ne_1S3"] = i;
      else if (level == "  1S2  ") mapLevels["Ne_1S2"] = i;
      else if (level == " 2P10  ") mapLevels["Ne_2P10"] = i;
      else if (level == "  2P9  ") mapLevels["Ne_2P9"] = i;
      else if (level == "  2P8  ") mapLevels["Ne_2P8"] = i;
      else if (level == "  2P7  ") mapLevels["Ne_2P7"] = i;
      else if (level == "  2P6  ") mapLevels["Ne_2P6"] = i;
      else if (level == "  2P5  ") mapLevels["Ne_2P5"] = i;
      else if (level == "  2P4  ") mapLevels["Ne_2P4"] = i;
      else if (level == "  2P3  ") mapLevels["Ne_2P3"] = i;
      else if (level == "  2P2  ") mapLevels["Ne_2P2"] = i;
      else if (level == "  2P1  ") mapLevels["Ne_2P1"] = i;
      else if (level == "  2S5  ") mapLevels["Ne_2S5"] = i;
      else if (level == "  2S4  ") mapLevels["Ne_2S4"] = i;
      else if (level == "  2S3  ") mapLevels["Ne_2S3"] = i;
      else if (level == "  2S2  ") mapLevels["Ne_2S2"] = i;
      else if (level == "  3D6  ") mapLevels["Ne_3D6"] = i;
      else if (level == "  3D5  ") mapLevels["Ne_3D5"] = i;
      else if (level == " 3D4!  ") mapLevels["Ne_3D4!"] = i;
      else if (level == "  3D4  ") mapLevels["Ne_3D4"] = i;
      else if (level == "  3D3  ") mapLevels["Ne_3D3"] = i;
      else if (level == "  3D2  ") mapLevels["Ne_3D2"] = i;
      else if (level == " 3D1!! ") mapLevels["Ne_3D1!!"] = i;
      else if (level == " 3D1!  ") mapLevels["Ne_3D1!"] = i;
      else if (level == "3S1!!!!") mapLevels["Ne_3S1!!!!"] = i;
      else if (level == "3S1!!! ") mapLevels["Ne_3S1!!!"] = i;
      else if (level == " 3S1!! ") mapLevels["Ne_3S1!!"] = i;
      else if (level == "  3S1! ") mapLevels["Ne_3S1!"] = i;
      else if (level == "SUM 3P1") mapLevels["Ne_3P10_3P6"] = i;
      else if (level == "SUM 3P5") mapLevels["Ne_3P5_3P2"] = i;
      else if (level == "  3P1  ") mapLevels["Ne_3P1"] = i;
      else if (level == "  3S4  ") mapLevels["Ne_3S4"] = i;
      else if (level == "  3S2  ") mapLevels["Ne_3S2"] = i;
      else if (level == "  4D5  ") mapLevels["Ne_4D5"] = i;
      else if (level == "  4D2  ") mapLevels["Ne_4D2"] = i;
      else if (level == " 4S1!  ") mapLevels["Ne_4S1!"] = i;
      else if (level == "  4S4  ") mapLevels["Ne_4S4"] = i;
      else if (level == "  5D5  ") mapLevels["Ne_5D5"] = i;
      else if (level == "  5D2  ") mapLevels["Ne_5D2"] = i;
      else if (level == "  4S2  ") mapLevels["Ne_4S2"] = i;
      else if (level == " 5S1!  ") mapLevels["Ne_5S1!"] = i;
      else if (level == "SUM S H") mapLevels["Ne_Sum_S_High"] = i;
      else if (level == "SUM D H") mapLevels["Ne_Sum_P_High"] = i;
      else {
        std::cerr << m_className << "::ComputeDeexcitationTable:\n";
        std::cerr << "    Unknown excitation level:\n";
        std::cerr << "      Ne " << level << "\n";
      }
      break;
    }
  }

  // Count the excitation levels.
  unsigned int nDeexcitations = 0;
  std::map<std::string, int> mapDxc;
  std::map<std::string, int>::iterator itMap;
  for (itMap = mapLevels.begin(); itMap != mapLevels.end(); itMap++) {
    std::string level = (*itMap).first;
    mapDxc[level] = nDeexcitations;
    m_iDeexcitation[(*itMap).second] = nDeexcitations;
    ++nDeexcitations;
  }

  // Conversion factor from oscillator strength to transition rate.
  const double f2A =
      2. * SpeedOfLight * FineStructureConstant / (3. * ElectronMass * HbarC);

  // Radiative de-excitation channels
  // Transition rates (unless indicated otherwise) are taken from:
  //     NIST Atomic Spectra Database
  // Transition rates for lines missing in the NIST database:
  //     O. Zatsarinny and K. Bartschat, J. Phys. B 39 (2006), 2145-2158
  // Oscillator strengths not included in the NIST database:
  //     J. Berkowitz, Atomic and Molecular Photoabsorption (2002)
  //     C.-M. Lee and K. T. Lu, Phys. Rev. A 8 (1973), 1241-1257
  for (itMap = mapLevels.begin(); itMap != mapLevels.end(); itMap++) {
    std::string level = (*itMap).first;
    deexcitation newDxc;
    newDxc.gas = int(m_csType[(*itMap).second] / nCsTypes);
    newDxc.level = (*itMap).second;
    newDxc.label = level;
    // Excitation energy
    newDxc.energy = m_energyLoss[(*itMap).second] * m_rgas[newDxc.gas];
    // Oscillator strength
    newDxc.osc = newDxc.cf = 0.;
    newDxc.sDoppler = newDxc.gPressure = newDxc.width = 0.;
    newDxc.nChannels = 0;
    if (level == "Ar_1S5" || level == "Ar_1S3") {
      // Metastables
    } else if (level == "Ar_1S4") {
      // Oscillator strength from NIST database
      newDxc.osc = 0.0609;
      // Berkowitz: f = 0.058
      int nc = 1;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 0.119;
      newDxc.final[0] = -1;
    } else if (level == "Ar_1S2") {
      // Oscillator strength from NIST database
      newDxc.osc = 0.25;
      // Berkowitz: 0.2214
      int nc = 1;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 0.51;
      newDxc.final[0] = -1;
    } else if (level == "Ar_2P10") {
      int nc = 4;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 0.0189;
      newDxc.final[0] = mapDxc["Ar_1S5"];
      newDxc.p[1] = 5.43e-3;
      newDxc.final[1] = mapDxc["Ar_1S4"];
      newDxc.p[2] = 9.8e-4;
      newDxc.final[2] = mapDxc["Ar_1S3"];
      newDxc.p[3] = 1.9e-4;
      newDxc.final[3] = mapDxc["Ar_1S2"];
    } else if (level == "Ar_2P9") {
      int nc = 1;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 0.0331;
      newDxc.final[0] = mapDxc["Ar_1S5"];
    } else if (level == "Ar_2P8") {
      int nc = 3;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 9.28e-3;
      newDxc.final[0] = mapDxc["Ar_1S5"];
      newDxc.p[1] = 0.0215;
      newDxc.final[1] = mapDxc["Ar_1S4"];
      newDxc.p[2] = 1.47e-3;
      newDxc.final[2] = mapDxc["Ar_1S2"];
    } else if (level == "Ar_2P7") {
      int nc = 4;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 5.18e-3;
      newDxc.final[0] = mapDxc["Ar_1S5"];
      newDxc.p[1] = 0.025;
      newDxc.final[1] = mapDxc["Ar_1S4"];
      newDxc.p[2] = 2.43e-3;
      newDxc.final[2] = mapDxc["Ar_1S3"];
      newDxc.p[3] = 1.06e-3;
      newDxc.final[3] = mapDxc["Ar_1S2"];
    } else if (level == "Ar_2P6") {
      int nc = 3;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 0.0245;
      newDxc.final[0] = mapDxc["Ar_1S5"];
      newDxc.p[1] = 4.9e-3;
      newDxc.final[1] = mapDxc["Ar_1S4"];
      newDxc.p[2] = 5.03e-3;
      newDxc.final[2] = mapDxc["Ar_1S2"];
    } else if (level == "Ar_2P5") {
      int nc = 1;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 0.0402;
      newDxc.final[0] = mapDxc["Ar_1S4"];
    } else if (level == "Ar_2P4") {
      int nc = 4;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 6.25e-4;
      newDxc.final[0] = mapDxc["Ar_1S5"];
      newDxc.p[1] = 2.2e-5;
      newDxc.final[1] = mapDxc["Ar_1S4"];
      newDxc.p[2] = 0.0186;
      newDxc.final[2] = mapDxc["Ar_1S3"];
      newDxc.p[3] = 0.0139;
      newDxc.final[3] = mapDxc["Ar_1S2"];
    } else if (level == "Ar_2P3") {
      int nc = 3;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 3.8e-3;
      newDxc.final[0] = mapDxc["Ar_1S5"];
      newDxc.p[1] = 8.47e-3;
      newDxc.final[1] = mapDxc["Ar_1S4"];
      newDxc.p[2] = 0.0223;
      newDxc.final[2] = mapDxc["Ar_1S2"];
    } else if (level == "Ar_2P2") {
      int nc = 4;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 6.39e-3;
      newDxc.final[0] = mapDxc["Ar_1S5"];
      newDxc.p[1] = 1.83e-3;
      newDxc.final[1] = mapDxc["Ar_1S4"];
      newDxc.p[2] = 0.0117;
      newDxc.final[2] = mapDxc["Ar_1S3"];
      newDxc.p[3] = 0.0153;
      newDxc.final[3] = mapDxc["Ar_1S2"];
    } else if (level == "Ar_2P1") {
      int nc = 2;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 2.36e-4;
      newDxc.final[0] = mapDxc["Ar_1S4"];
      newDxc.p[1] = 0.0445;
      newDxc.final[1] = mapDxc["Ar_1S2"];
    } else if (level == "Ar_3D6") {
      int nc = 4;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional line (2P7) from Bartschat
      newDxc.p[0] = 8.1e-3;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 7.73e-4;
      newDxc.final[1] = mapDxc["Ar_2P7"];
      newDxc.p[2] = 1.2e-4;
      newDxc.final[2] = mapDxc["Ar_2P4"];
      newDxc.p[3] = 3.6e-4;
      newDxc.final[3] = mapDxc["Ar_2P2"];
    } else if (level == "Ar_3D5") {
      // Oscillator strength from Berkowitz
      newDxc.osc = 0.0011;
      int nc = 10;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional lines (2P7, 2P6, 2P5, 2P1) from Bartschat
      newDxc.p[0] = 7.4e-3;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 3.9e-5;
      newDxc.final[1] = mapDxc["Ar_2P8"];
      newDxc.p[2] = 3.09e-4;
      newDxc.final[2] = mapDxc["Ar_2P7"];
      newDxc.p[3] = 1.37e-3;
      newDxc.final[3] = mapDxc["Ar_2P6"];
      newDxc.p[4] = 5.75e-4;
      newDxc.final[4] = mapDxc["Ar_2P5"];
      newDxc.p[5] = 3.2e-5;
      newDxc.final[5] = mapDxc["Ar_2P4"];
      newDxc.p[6] = 1.4e-4;
      newDxc.final[6] = mapDxc["Ar_2P3"];
      newDxc.p[7] = 1.7e-4;
      newDxc.final[7] = mapDxc["Ar_2P2"];
      newDxc.p[8] = 2.49e-6;
      newDxc.final[8] = mapDxc["Ar_2P1"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[9] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[9] = -1;
    } else if (level == "Ar_3D3") {
      int nc = 8;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional lines (2P9, 2P4) from Bartschat
      newDxc.p[0] = 4.9e-3;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 9.82e-5;
      newDxc.final[1] = mapDxc["Ar_2P9"];
      newDxc.p[2] = 1.2e-4;
      newDxc.final[2] = mapDxc["Ar_2P8"];
      newDxc.p[3] = 2.6e-4;
      newDxc.final[3] = mapDxc["Ar_2P7"];
      newDxc.p[4] = 2.5e-3;
      newDxc.final[4] = mapDxc["Ar_2P6"];
      newDxc.p[5] = 9.41e-5;
      newDxc.final[5] = mapDxc["Ar_2P4"];
      newDxc.p[6] = 3.9e-4;
      newDxc.final[6] = mapDxc["Ar_2P3"];
      newDxc.p[7] = 1.1e-4;
      newDxc.final[7] = mapDxc["Ar_2P2"];
    } else if (level == "Ar_3D4!") {
      int nc = 1;
      newDxc.nChannels = nc;
      // Transition probability for 2P9 transition from Bartschat
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 0.01593;
      newDxc.final[0] = mapDxc["Ar_2P9"];
    } else if (level == "Ar_3D4") {
      int nc = 4;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional lines (2P9, 2P3) from Bartschat
      newDxc.p[0] = 2.29e-3;
      newDxc.final[0] = mapDxc["Ar_2P9"];
      newDxc.p[1] = 0.011;
      newDxc.final[1] = mapDxc["Ar_2P8"];
      newDxc.p[2] = 8.8e-5;
      newDxc.final[2] = mapDxc["Ar_2P6"];
      newDxc.p[3] = 2.53e-6;
      newDxc.final[3] = mapDxc["Ar_2P3"];
    } else if (level == "Ar_3D1!!") {
      int nc = 8;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional lines (2P10, 2P6, 2P4 - 2P2) from Bartschat
      newDxc.p[0] = 5.85e-6;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 1.2e-4;
      newDxc.final[1] = mapDxc["Ar_2P9"];
      newDxc.p[2] = 5.7e-3;
      newDxc.final[2] = mapDxc["Ar_2P8"];
      newDxc.p[3] = 7.3e-3;
      newDxc.final[3] = mapDxc["Ar_2P7"];
      newDxc.p[4] = 2.e-4;
      newDxc.final[4] = mapDxc["Ar_2P6"];
      newDxc.p[5] = 1.54e-6;
      newDxc.final[5] = mapDxc["Ar_2P4"];
      newDxc.p[6] = 2.08e-5;
      newDxc.final[6] = mapDxc["Ar_2P3"];
      newDxc.p[7] = 6.75e-7;
      newDxc.final[7] = mapDxc["Ar_2P2"];
    } else if (level == "Ar_2S5") {
      int nc = 8;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 4.9e-3;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 0.011;
      newDxc.final[1] = mapDxc["Ar_2P9"];
      newDxc.p[2] = 1.1e-3;
      newDxc.final[2] = mapDxc["Ar_2P8"];
      newDxc.p[3] = 4.6e-4;
      newDxc.final[3] = mapDxc["Ar_2P7"];
      newDxc.p[4] = 3.3e-3;
      newDxc.final[4] = mapDxc["Ar_2P6"];
      newDxc.p[5] = 5.9e-5;
      newDxc.final[5] = mapDxc["Ar_2P4"];
      newDxc.p[6] = 1.2e-4;
      newDxc.final[6] = mapDxc["Ar_2P3"];
      newDxc.p[7] = 3.1e-4;
      newDxc.final[7] = mapDxc["Ar_2P2"];
    } else if (level == "Ar_2S4") {
      // Oscillator strength from NIST database
      newDxc.osc = 0.027;
      // Berkowitz: f = 0.026;
      int nc = 10;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 0.077;
      newDxc.final[0] = -1;
      newDxc.p[1] = 2.44e-3;
      newDxc.final[1] = mapDxc["Ar_2P10"];
      newDxc.p[2] = 8.9e-3;
      newDxc.final[2] = mapDxc["Ar_2P8"];
      newDxc.p[3] = 4.6e-3;
      newDxc.final[3] = mapDxc["Ar_2P7"];
      newDxc.p[4] = 2.7e-3;
      newDxc.final[4] = mapDxc["Ar_2P6"];
      newDxc.p[5] = 1.3e-3;
      newDxc.final[5] = mapDxc["Ar_2P5"];
      newDxc.p[6] = 4.5e-4;
      newDxc.final[6] = mapDxc["Ar_2P4"];
      newDxc.p[7] = 2.9e-5;
      newDxc.final[7] = mapDxc["Ar_2P3"];
      newDxc.p[8] = 3.e-5;
      newDxc.final[8] = mapDxc["Ar_2P2"];
      newDxc.p[9] = 1.6e-4;
      newDxc.final[9] = mapDxc["Ar_2P1"];
    } else if (level == "Ar_3D1!") {
      int nc = 4;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional line (2P6) from Bartschat
      newDxc.p[0] = 3.1e-3;
      newDxc.final[0] = mapDxc["Ar_2P9"];
      newDxc.p[1] = 2.e-3;
      newDxc.final[1] = mapDxc["Ar_2P8"];
      newDxc.p[2] = 0.015;
      newDxc.final[2] = mapDxc["Ar_2P6"];
      newDxc.p[3] = 9.8e-6;
      newDxc.final[3] = mapDxc["Ar_2P3"];
    } else if (level == "Ar_3D2") {
      // Oscillator strength from NIST database
      newDxc.osc = 0.0932;
      // Berkowitz: f = 0.09
      int nc = 10;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional lines (2P10, 2P6, 2P4-2P1) from Bartschat
      newDxc.p[0] = 0.27;
      newDxc.final[0] = -1;
      newDxc.p[1] = 1.35e-5;
      newDxc.final[1] = mapDxc["Ar_2P10"];
      newDxc.p[2] = 9.52e-4;
      newDxc.final[2] = mapDxc["Ar_2P8"];
      newDxc.p[3] = 0.011;
      newDxc.final[3] = mapDxc["Ar_2P7"];
      newDxc.p[4] = 4.01e-5;
      newDxc.final[4] = mapDxc["Ar_2P6"];
      newDxc.p[5] = 4.3e-3;
      newDxc.final[5] = mapDxc["Ar_2P5"];
      newDxc.p[6] = 8.96e-4;
      newDxc.final[6] = mapDxc["Ar_2P4"];
      newDxc.p[7] = 4.45e-5;
      newDxc.final[7] = mapDxc["Ar_2P3"];
      newDxc.p[8] = 5.87e-5;
      newDxc.final[8] = mapDxc["Ar_2P2"];
      newDxc.p[9] = 8.77e-4;
      newDxc.final[9] = mapDxc["Ar_2P1"];
    } else if (level == "Ar_3S1!!!!") {
      int nc = 8;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional lines (2P10, 2P9, 2P7, 2P6, 2P2) from Bartschat
      newDxc.p[0] = 7.51e-6;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 4.3e-5;
      newDxc.final[1] = mapDxc["Ar_2P9"];
      newDxc.p[2] = 8.3e-4;
      newDxc.final[2] = mapDxc["Ar_2P8"];
      newDxc.p[3] = 5.01e-5;
      newDxc.final[3] = mapDxc["Ar_2P7"];
      newDxc.p[4] = 2.09e-4;
      newDxc.final[4] = mapDxc["Ar_2P6"];
      newDxc.p[5] = 0.013;
      newDxc.final[5] = mapDxc["Ar_2P4"];
      newDxc.p[6] = 2.2e-3;
      newDxc.final[6] = mapDxc["Ar_2P3"];
      newDxc.p[7] = 3.35e-6;
      newDxc.final[7] = mapDxc["Ar_2P2"];
    } else if (level == "Ar_3S1!!") {
      int nc = 8;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional lines (2P10 - 2P8, 2P4, 2P3)
      newDxc.p[0] = 1.89e-4;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 1.52e-4;
      newDxc.final[1] = mapDxc["Ar_2P9"];
      newDxc.p[2] = 7.21e-4;
      newDxc.final[2] = mapDxc["Ar_2P8"];
      newDxc.p[3] = 3.69e-4;
      newDxc.final[3] = mapDxc["Ar_2P7"];
      newDxc.p[4] = 3.76e-3;
      newDxc.final[4] = mapDxc["Ar_2P6"];
      newDxc.p[5] = 1.72e-4;
      newDxc.final[5] = mapDxc["Ar_2P4"];
      newDxc.p[6] = 5.8e-4;
      newDxc.final[6] = mapDxc["Ar_2P3"];
      newDxc.p[7] = 6.2e-3;
      newDxc.final[7] = mapDxc["Ar_2P2"];
    } else if (level == "Ar_3S1!!!") {
      int nc = 4;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional lines (2P9, 2P8, 2P6) from Bartschat
      newDxc.p[0] = 7.36e-4;
      newDxc.final[0] = mapDxc["Ar_2P9"];
      newDxc.p[1] = 4.2e-5;
      newDxc.final[1] = mapDxc["Ar_2P8"];
      newDxc.p[2] = 9.3e-5;
      newDxc.final[2] = mapDxc["Ar_2P6"];
      newDxc.p[3] = 0.015;
      newDxc.final[3] = mapDxc["Ar_2P3"];
    } else if (level == "Ar_2S3") {
      int nc = 4;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 3.26e-3;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 2.22e-3;
      newDxc.final[1] = mapDxc["Ar_2P7"];
      newDxc.p[2] = 0.01;
      newDxc.final[2] = mapDxc["Ar_2P4"];
      newDxc.p[3] = 5.1e-3;
      newDxc.final[3] = mapDxc["Ar_2P2"];
    } else if (level == "Ar_2S2") {
      // Oscillator strength from NIST database
      newDxc.osc = 0.0119;
      // Berkowitz: f = 0.012;
      int nc = 10;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 0.035;
      newDxc.final[0] = -1;
      newDxc.p[1] = 1.76e-3;
      newDxc.final[1] = mapDxc["Ar_2P10"];
      newDxc.p[2] = 2.1e-4;
      newDxc.final[2] = mapDxc["Ar_2P8"];
      newDxc.p[3] = 2.8e-4;
      newDxc.final[3] = mapDxc["Ar_2P7"];
      newDxc.p[4] = 1.39e-3;
      newDxc.final[4] = mapDxc["Ar_2P6"];
      newDxc.p[5] = 3.8e-4;
      newDxc.final[5] = mapDxc["Ar_2P5"];
      newDxc.p[6] = 2.0e-3;
      newDxc.final[6] = mapDxc["Ar_2P4"];
      newDxc.p[7] = 8.9e-3;
      newDxc.final[7] = mapDxc["Ar_2P3"];
      newDxc.p[8] = 3.4e-3;
      newDxc.final[8] = mapDxc["Ar_2P2"];
      newDxc.p[9] = 1.9e-3;
      newDxc.final[9] = mapDxc["Ar_2P1"];
    } else if (level == "Ar_3S1!") {
      // Oscillator strength from NIST database
      newDxc.osc = 0.106;
      // Berkowitz: f = 0.106
      int nc = 10;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional lines (2P10, 2P8, 2P7, 2P3) from Bartschat
      newDxc.p[0] = 0.313;
      newDxc.final[0] = -1;
      newDxc.p[1] = 2.05e-5;
      newDxc.final[1] = mapDxc["Ar_2P10"];
      newDxc.p[2] = 8.33e-5;
      newDxc.final[2] = mapDxc["Ar_2P8"];
      newDxc.p[3] = 3.9e-4;
      newDxc.final[3] = mapDxc["Ar_2P7"];
      newDxc.p[4] = 3.96e-4;
      newDxc.final[4] = mapDxc["Ar_2P6"];
      newDxc.p[5] = 4.2e-4;
      newDxc.final[5] = mapDxc["Ar_2P5"];
      newDxc.p[6] = 4.5e-3;
      newDxc.final[6] = mapDxc["Ar_2P4"];
      newDxc.p[7] = 4.84e-5;
      newDxc.final[7] = mapDxc["Ar_2P3"];
      newDxc.p[8] = 7.1e-3;
      newDxc.final[8] = mapDxc["Ar_2P2"];
      newDxc.p[9] = 5.2e-3;
      newDxc.final[9] = mapDxc["Ar_2P1"];
    } else if (level == "Ar_4D5") {
      // Oscillator strength from Berkowitz
      newDxc.osc = 0.0019;
      int nc = 7;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 2.78e-3;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 2.8e-4;
      newDxc.final[1] = mapDxc["Ar_2P8"];
      newDxc.p[2] = 8.6e-4;
      newDxc.final[2] = mapDxc["Ar_2P6"];
      newDxc.p[3] = 9.2e-4;
      newDxc.final[3] = mapDxc["Ar_2P5"];
      newDxc.p[4] = 4.6e-4;
      newDxc.final[4] = mapDxc["Ar_2P3"];
      newDxc.p[5] = 1.6e-4;
      newDxc.final[5] = mapDxc["Ar_2P2"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[6] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[6] = -1;
    } else if (level == "Ar_3S4") {
      // Oscillator strength from Berkowitz
      newDxc.osc = 0.0144;
      int nc = 10;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 4.21e-4;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 2.e-3;
      newDxc.final[1] = mapDxc["Ar_2P8"];
      newDxc.p[2] = 1.7e-3;
      newDxc.final[2] = mapDxc["Ar_2P7"];
      newDxc.p[3] = 7.2e-4;
      newDxc.final[3] = mapDxc["Ar_2P6"];
      newDxc.p[4] = 3.5e-4;
      newDxc.final[4] = mapDxc["Ar_2P5"];
      newDxc.p[5] = 1.2e-4;
      newDxc.final[5] = mapDxc["Ar_2P4"];
      newDxc.p[6] = 4.2e-6;
      newDxc.final[6] = mapDxc["Ar_2P3"];
      newDxc.p[7] = 3.3e-5;
      newDxc.final[7] = mapDxc["Ar_2P2"];
      newDxc.p[8] = 9.7e-5;
      newDxc.final[8] = mapDxc["Ar_2P1"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[9] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[9] = -1;
    } else if (level == "Ar_4D2") {
      // Oscillator strength from Berkowitz
      newDxc.osc = 0.048;
      int nc = 2;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 1.7e-4;
      newDxc.final[0] = mapDxc["Ar_2P7"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[1] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[1] = -1;
    } else if (level == "Ar_4S1!") {
      // Oscillator strength from Berkowitz
      newDxc.osc = 0.0209;
      int nc = 7;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 1.05e-3;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 3.1e-5;
      newDxc.final[1] = mapDxc["Ar_2P8"];
      newDxc.p[2] = 2.5e-5;
      newDxc.final[2] = mapDxc["Ar_2P7"];
      newDxc.p[3] = 4.0e-4;
      newDxc.final[3] = mapDxc["Ar_2P6"];
      newDxc.p[4] = 5.8e-5;
      newDxc.final[4] = mapDxc["Ar_2P5"];
      newDxc.p[5] = 1.2e-4;
      newDxc.final[5] = mapDxc["Ar_2P3"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[6] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[6] = -1;
    } else if (level == "Ar_3S2") {
      // Oscillator strength from Berkowitz
      newDxc.osc = 0.0221;
      int nc = 10;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 2.85e-4;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 5.1e-5;
      newDxc.final[1] = mapDxc["Ar_2P8"];
      newDxc.p[2] = 5.3e-5;
      newDxc.final[2] = mapDxc["Ar_2P7"];
      newDxc.p[3] = 1.6e-4;
      newDxc.final[3] = mapDxc["Ar_2P6"];
      newDxc.p[4] = 1.5e-4;
      newDxc.final[4] = mapDxc["Ar_2P5"];
      newDxc.p[5] = 6.0e-4;
      newDxc.final[5] = mapDxc["Ar_2P4"];
      newDxc.p[6] = 2.48e-3;
      newDxc.final[6] = mapDxc["Ar_2P3"];
      newDxc.p[7] = 9.6e-4;
      newDxc.final[7] = mapDxc["Ar_2P2"];
      newDxc.p[8] = 3.59e-4;
      newDxc.final[8] = mapDxc["Ar_2P1"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[9] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[9] = -1;
    } else if (level == "Ar_5D5") {
      // Oscillator strength from Berkowitz
      newDxc.osc = 0.0041;
      int nc = 9;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 2.2e-3;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 1.1e-4;
      newDxc.final[1] = mapDxc["Ar_2P8"];
      newDxc.p[2] = 7.6e-5;
      newDxc.final[2] = mapDxc["Ar_2P7"];
      newDxc.p[3] = 4.2e-4;
      newDxc.final[3] = mapDxc["Ar_2P6"];
      newDxc.p[4] = 2.4e-4;
      newDxc.final[4] = mapDxc["Ar_2P5"];
      newDxc.p[5] = 2.1e-4;
      newDxc.final[5] = mapDxc["Ar_2P4"];
      newDxc.p[6] = 2.4e-4;
      newDxc.final[6] = mapDxc["Ar_2P3"];
      newDxc.p[7] = 1.2e-4;
      newDxc.final[7] = mapDxc["Ar_2P2"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[8] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[8] = -1;
    } else if (level == "Ar_4S4") {
      // Oscillator strength from Berkowitz
      newDxc.osc = 0.0139;
      int nc = 7;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 1.9e-4;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 1.1e-3;
      newDxc.final[1] = mapDxc["Ar_2P8"];
      newDxc.p[2] = 5.2e-4;
      newDxc.final[2] = mapDxc["Ar_2P7"];
      newDxc.p[3] = 5.1e-4;
      newDxc.final[3] = mapDxc["Ar_2P6"];
      newDxc.p[4] = 9.4e-5;
      newDxc.final[4] = mapDxc["Ar_2P5"];
      newDxc.p[5] = 5.4e-5;
      newDxc.final[5] = mapDxc["Ar_2P4"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[6] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[6] = -1;
    } else if (level == "Ar_5D2") {
      // Oscillator strength from Berkowitz
      newDxc.osc = 0.0426;
      int nc = 5;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 5.9e-5;
      newDxc.final[0] = mapDxc["Ar_2P8"];
      newDxc.p[1] = 9.0e-6;
      newDxc.final[1] = mapDxc["Ar_2P7"];
      newDxc.p[2] = 1.5e-4;
      newDxc.final[2] = mapDxc["Ar_2P5"];
      newDxc.p[3] = 3.1e-5;
      newDxc.final[3] = mapDxc["Ar_2P2"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[4] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[4] = -1;
    } else if (level == "Ar_6D5") {
      // Oscillator strength from Lee and Lu
      newDxc.osc = 0.00075;
      // Berkowitz estimates f = 0.0062 for the sum of
      // all "weak" nd levels with n = 6 and higher.
      int nc = 7;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 1.9e-3;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 4.2e-4;
      newDxc.final[1] = mapDxc["Ar_2P6"];
      newDxc.p[2] = 3.e-4;
      newDxc.final[2] = mapDxc["Ar_2P5"];
      newDxc.p[3] = 5.1e-5;
      newDxc.final[3] = mapDxc["Ar_2P4"];
      newDxc.p[4] = 6.6e-5;
      newDxc.final[4] = mapDxc["Ar_2P3"];
      newDxc.p[5] = 1.21e-4;
      newDxc.final[5] = mapDxc["Ar_2P1"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[6] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[6] = -1;
    } else if (level == "Ar_5S1!") {
      // Oscillator strength from Lee and Lu
      newDxc.osc = 0.00051;
      // Berkowitz estimates f = 0.0562 for the sum
      // of all nd' levels with n = 5 and higher.
      int nc = 2;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 7.7e-5;
      newDxc.final[0] = mapDxc["Ar_2P5"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[1] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[1] = -1;
    } else if (level == "Ar_4S2") {
      // Oscillator strength from Lee and Lu
      newDxc.osc = 0.00074;
      // Berkowitz estimates f = 0.0069 for the sum over all
      // ns' levels with n = 7 and higher.
      int nc = 8;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 4.5e-4;
      newDxc.final[0] = mapDxc["Ar_2P10"];
      newDxc.p[1] = 2.e-4;
      newDxc.final[1] = mapDxc["Ar_2P8"];
      newDxc.p[2] = 2.1e-4;
      newDxc.final[2] = mapDxc["Ar_2P7"];
      newDxc.p[3] = 1.2e-4;
      newDxc.final[3] = mapDxc["Ar_2P5"];
      newDxc.p[4] = 1.8e-4;
      newDxc.final[4] = mapDxc["Ar_2P4"];
      newDxc.p[5] = 9.e-4;
      newDxc.final[5] = mapDxc["Ar_2P3"];
      newDxc.p[6] = 3.3e-4;
      newDxc.final[6] = mapDxc["Ar_2P2"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[7] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[7] = -1;
    } else if (level == "Ar_5S4") {
      // Oscillator strength from Lee and Lu
      newDxc.osc = 0.0130;
      // Berkowitz estimates f = 0.0211 for the sum of all
      // ns levels with n = 8 and higher.
      newDxc.osc = 0.0211;
      int nc = 6;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      newDxc.p[0] = 3.6e-4;
      newDxc.final[0] = mapDxc["Ar_2P8"];
      newDxc.p[1] = 1.2e-4;
      newDxc.final[1] = mapDxc["Ar_2P6"];
      newDxc.p[2] = 1.5e-4;
      newDxc.final[2] = mapDxc["Ar_2P4"];
      newDxc.p[3] = 1.4e-4;
      newDxc.final[3] = mapDxc["Ar_2P3"];
      newDxc.p[4] = 7.5e-5;
      newDxc.final[4] = mapDxc["Ar_2P2"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[5] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[5] = -1;
    } else if (level == "Ar_6D2") {
      // Oscillator strength from Lee and Lu
      newDxc.osc = 0.0290;
      // Berkowitz estimates f = 0.0574 for the sum of all
      // "strong" nd levels with n = 6 and higher.
      newDxc.osc = 0.0574;
      int nc = 2;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeRad);
      // Additional line: 2P7
      newDxc.p[0] = 3.33e-3;
      newDxc.final[0] = mapDxc["Ar_2P7"];
      // Transition probability to ground state calculated from osc. strength
      newDxc.p[1] = f2A * pow(newDxc.energy, 2) * newDxc.osc;
      newDxc.final[1] = -1;
    } else if (level == "Ar_Higher") {
      newDxc.osc = 0.;
      // This (artificial) level represents the sum of higher J = 1 states.
      // The deeexcitation cascade is simulated by allocating it
      // with equal probability to one of the five nearest levels below.
      int nc = 5;
      newDxc.nChannels = nc;
      newDxc.p.resize(nc);
      newDxc.final.resize(nc);
      newDxc.type.resize(nc, DxcTypeCollNonIon);
      newDxc.p[0] = 100.;
      newDxc.final[0] = mapDxc["Ar_6D5"];
      newDxc.p[1] = 100.;
      newDxc.final[1] = mapDxc["Ar_5S1!"];
      newDxc.p[2] = 100.;
      newDxc.final[2] = mapDxc["Ar_4S2"];
      newDxc.p[3] = 100.;
      newDxc.final[3] = mapDxc["Ar_5S4"];
      newDxc.p[4] = 100.;
      newDxc.final[4] = mapDxc["Ar_6D2"];
    } else {
      std::cerr << m_className << "::ComputeDeexcitationTable:\n";
      std::cerr << "    Missing de-excitation data for level " << level
                << ".\n";
      std::cerr << "    Program bug!\n";
      return;
    }
    m_deexcitations.push_back(newDxc);
  }

  if (m_debug || verbose) {
    std::cout << m_className << "::ComputeDeexcitationTable:\n";
    std::cout << "    Found " << m_deexcitations.size() << " levels "
              << "with available radiative de-excitation data.\n";
  }

  // Collisional de-excitation channels
  if (withAr) {
    // Add the Ar dimer ground state.
    deexcitation dimer;
    dimer.label = "Ar_Dimer";
    dimer.level = -1;
    dimer.gas = iAr;
    dimer.energy = 14.71;
    dimer.osc = dimer.cf = 0.;
    dimer.sDoppler = dimer.gPressure = dimer.width = 0.;
    dimer.nChannels = 0;
    mapDxc["Ar_Dimer"] = m_deexcitations.size();
    m_deexcitations.push_back(dimer);
    ++nDeexcitations;
    // Add an Ar excimer level.
    deexcitation excimer;
    excimer.label = "Ar_Excimer";
    excimer.level = -1;
    excimer.gas = iAr;
    excimer.energy = 14.71;
    excimer.osc = excimer.cf = 0.;
    excimer.sDoppler = excimer.gPressure = excimer.width = 0.;
    excimer.nChannels = 0;
    mapDxc["Ar_Excimer"] = m_deexcitations.size();
    m_deexcitations.push_back(excimer);
    ++nDeexcitations;
    const double nAr = GetNumberDensity() * cAr;
    for (unsigned int j = 0; j < nDeexcitations; ++j) {
      const std::string level = m_deexcitations[j].label;
      if (level == "Ar_1S5") {
        // Two-body and three-body collision rate constants
        // Three-body collisions lead to excimer formation.
        // Two-body collisions give rise to collisional mixing.
        const bool useTachibanaData = false;
        const bool useKoltsSetserData = true;
        const bool useCollMixing = true;
        if (useTachibanaData) {
          // K. Tachibana, Phys. Rev. A 34 (1986), 1007-1015
          const double k2b = 2.3e-24;
          const double k3b = 1.4e-41;
          m_deexcitations[j].p.push_back(k3b * nAr * nAr);
          m_deexcitations[j].final.push_back(mapDxc["Ar_Excimer"]);
          m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
          m_deexcitations[j].nChannels += 1;
          if (useCollMixing) {
            m_deexcitations[j].p.push_back(k2b * nAr);
            m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
            m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
            m_deexcitations[j].nChannels += 1;
          }
        } else if (useKoltsSetserData) {
          // Kolts and Setser, J. Chem. Phys. 68 (1978), 4848-4859
          const double k2b = 2.1e-24;
          const double k3b = 1.1e-41;
          m_deexcitations[j].p.push_back(k3b * nAr * nAr);
          m_deexcitations[j].final.push_back(mapDxc["Ar_Excimer"]);
          m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
          m_deexcitations[j].nChannels += 1;
          if (useCollMixing) {
            m_deexcitations[j].p.push_back(k2b * nAr);
            m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
            m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
            m_deexcitations[j].nChannels += 1;
          }
        }
      }
      if (level == "Ar_1S3") {
        // Two-body and three-body collision rate constants
        const bool useTachibanaData = false;
        const bool useKoltsSetserData = true;
        const bool useCollMixing = true;
        if (useTachibanaData) {
          // K. Tachibana, Phys. Rev. A 34 (1986), 1007-1015
          const double k2b = 4.3e-24;
          const double k3b = 1.5e-41;
          m_deexcitations[j].p.push_back(k3b * nAr * nAr);
          m_deexcitations[j].final.push_back(mapDxc["Ar_Excimer"]);
          m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
          m_deexcitations[j].nChannels += 1;
          if (useCollMixing) {
            m_deexcitations[j].p.push_back(k2b * nAr);
            m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
            m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
            m_deexcitations[j].nChannels += 1;
          }
        } else if (useKoltsSetserData) {
          // Kolts and Setser, J. Chem. Phys. 68 (1978), 4848-4859
          const double k2b = 5.3e-24;
          const double k3b = 0.83e-41;
          m_deexcitations[j].p.push_back(k3b * nAr * nAr);
          m_deexcitations[j].final.push_back(mapDxc["Ar_Excimer"]);
          m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
          m_deexcitations[j].nChannels += 1;
          if (useCollMixing) {
            m_deexcitations[j].p.push_back(k2b * nAr);
            m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
            m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
            m_deexcitations[j].nChannels += 1;
          }
        }
      }
      if (level == "Ar_2P1") {
        // Transfer to 4s states
        // Inoue, Setser, and Sadeghi, J. Chem. Phys. 75 (1982), 977-983
        // const double k4s = 2.9e-20;
        // Sadeghi et al. J. Chem. Phys. 115 (2001), 3144-3154
        const double k4s = 1.6e-20;
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S3"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S2"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 4;
      }
      if (level == "Ar_2P2") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        const double k23 = 0.5e-21;
        m_deexcitations[j].p.push_back(k23 * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P3"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
        // Transfer to 4s states
        // Inoue, Setser, and Sadeghi, J. Chem. Phys. 75 (1982), 977-983
        // const double k4s = 3.8e-20;
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        const double k4s = 5.3e-20;
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S3"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S2"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 4;
      }
      if (level == "Ar_2P3") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        const double k34 = 27.5e-21;
        const double k35 = 0.3e-21;
        const double k36 = 44.0e-21;
        const double k37 = 1.4e-21;
        const double k38 = 1.9e-21;
        const double k39 = 0.8e-21;
        m_deexcitations[j].p.push_back(k34 * nAr);
        m_deexcitations[j].p.push_back(k35 * nAr);
        m_deexcitations[j].p.push_back(k36 * nAr);
        m_deexcitations[j].p.push_back(k37 * nAr);
        m_deexcitations[j].p.push_back(k38 * nAr);
        m_deexcitations[j].p.push_back(k39 * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P6"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P7"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P8"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P9"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 6;
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        const double k4s = 4.7e-20;
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S3"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S2"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 4;
      }
      if (level == "Ar_2P4") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        const double k43 = 23.0e-21;
        const double k45 = 0.7e-21;
        const double k46 = 4.8e-21;
        const double k47 = 3.2e-21;
        const double k48 = 1.4e-21;
        const double k49 = 3.3e-21;
        m_deexcitations[j].p.push_back(k43 * nAr);
        m_deexcitations[j].p.push_back(k45 * nAr);
        m_deexcitations[j].p.push_back(k46 * nAr);
        m_deexcitations[j].p.push_back(k47 * nAr);
        m_deexcitations[j].p.push_back(k48 * nAr);
        m_deexcitations[j].p.push_back(k49 * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P3"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P6"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P7"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P8"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P9"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 6;
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        const double k4s = 3.9e-20;
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S3"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S2"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 4;
      }
      if (level == "Ar_2P5") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        const double k54 = 1.7e-21;
        const double k56 = 11.3e-21;
        const double k58 = 9.5e-21;
        m_deexcitations[j].p.push_back(k54 * nAr);
        m_deexcitations[j].p.push_back(k56 * nAr);
        m_deexcitations[j].p.push_back(k58 * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P6"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P8"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 3;
      }
      if (level == "Ar_2P6") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        const double k67 = 4.1e-21;
        const double k68 = 6.0e-21;
        const double k69 = 1.0e-21;
        m_deexcitations[j].p.push_back(k67 * nAr);
        m_deexcitations[j].p.push_back(k68 * nAr);
        m_deexcitations[j].p.push_back(k69 * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P7"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P8"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P9"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 3;
      }
      if (level == "Ar_2P7") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        const double k76 = 2.5e-21;
        const double k78 = 14.3e-21;
        const double k79 = 23.3e-21;
        m_deexcitations[j].p.push_back(k76 * nAr);
        m_deexcitations[j].p.push_back(k78 * nAr);
        m_deexcitations[j].p.push_back(k79 * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P6"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P8"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P9"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 3;
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        const double k4s = 5.5e-20;
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S3"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S2"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 4;
      }
      if (level == "Ar_2P8") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        const double k86 = 0.3e-21;
        const double k87 = 0.8e-21;
        const double k89 = 18.2e-21;
        const double k810 = 1.0e-21;
        m_deexcitations[j].p.push_back(k86 * nAr);
        m_deexcitations[j].p.push_back(k87 * nAr);
        m_deexcitations[j].p.push_back(k89 * nAr);
        m_deexcitations[j].p.push_back(k810 * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P6"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P7"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P9"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P10"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 4;
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        const double k4s = 3.e-20;
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S3"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S2"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 4;
      }
      if (level == "Ar_2P9") {
        // Collisional population transfer within 4p levels
        // T. D. Nguyen and N. Sadeghi, Phys. Rev. 18 (1978), 1388-1395
        const double k98 = 6.8e-21;
        const double k910 = 5.1e-21;
        m_deexcitations[j].p.push_back(k98 * nAr);
        m_deexcitations[j].p.push_back(k910 * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P8"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P10"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        const double k4s = 3.5e-20;
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S3"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S2"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 4;
      }
      if (level == "Ar_2P10") {
        // Transfer to 4s states
        // Chang and Setser, J. Chem. Phys. 69 (1978), 3885-3897
        const double k4s = 2.0e-20;
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].p.push_back(0.25 * k4s * nAr);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S3"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_1S2"]);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 4;
      }
      if (level == "Ar_3D6" || level == "Ar_3D5" || level == "Ar_3D3" ||
          level == "Ar_3D4!" || level == "Ar_3D4" || level == "Ar_3D1!!" ||
          level == "Ar_3D1!" || level == "Ar_3D2" || level == "Ar_3S1!!!!" ||
          level == "Ar_3S1!!" || level == "Ar_3S1!!!" || level == "Ar_3S1!" ||
          level == "Ar_2S5" || level == "Ar_2S4" || level == "Ar_2S3" ||
          level == "Ar_2S2") {
        // 3d and 5s levels
        // Transfer to 4p levels
        const double k4p = fit3d4p * 1.e-20;
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P10"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P9"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P8"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P7"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P6"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P3"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P2"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P1"]);
        for (int k = 10; k--;) {
          m_deexcitations[j].p.push_back(0.1 * k4p * nAr);
          m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        }
        m_deexcitations[j].nChannels += 10;
      }
      if (level == "Ar_4D5" || level == "Ar_3S4" || level == "Ar_4D2" ||
          level == "Ar_4S1!" || level == "Ar_3S2" || level == "Ar_5D5" ||
          level == "Ar_4S4" || level == "Ar_5D2" || level == "Ar_6D5" ||
          level == "Ar_5S1!" || level == "Ar_4S2" || level == "Ar_5S4" ||
          level == "Ar_6D2") {
        // Transfer to 4p levels
        const double k4p = fitHigh4p * 1.e-20;
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P10"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P9"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P8"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P7"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P6"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P5"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P4"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P3"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P2"]);
        m_deexcitations[j].final.push_back(mapDxc["Ar_2P1"]);
        for (int k = 10; k--;) {
          m_deexcitations[j].p.push_back(0.1 * k4p * nAr);
          m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        }
        m_deexcitations[j].nChannels += 10;
        // Hornbeck-Molnar ionisation
        // P. Becker and F. Lampe, J. Chem. Phys. 42 (1965), 3857-3863
        // A. Bogaerts and R. Gijbels, Phys. Rev. A 52 (1995), 3743-3751
        // This value seems high, to be checked!
        const double kHM = 2.e-18;
        const bool useHornbeckMolnar = true;
        if (useHornbeckMolnar) {
          m_deexcitations[j].p.push_back(kHM * nAr);
          m_deexcitations[j].final.push_back(mapDxc["Ar_Dimer"]);
          m_deexcitations[j].type.push_back(DxcTypeCollIon);
          m_deexcitations[j].nChannels += 1;
        }
      }
    }
  }

  // Collisional deexcitation by quenching gases.
  bool withCO2 = false;
  double cCO2 = 0.;
  int iCO2 = 0;
  bool withCH4 = false;
  double cCH4 = 0.;
  int iCH4 = 0;
  bool withC2H6 = false;
  double cC2H6 = 0.;
  int iC2H6 = 0;
  bool withIso = false;
  double cIso = 0.;
  int iIso = 0;
  bool withC2H2 = false;
  double cC2H2 = 0.;
  int iC2H2 = 0;
  bool withCF4 = false;
  double cCF4 = 0.;
  int iCF4 = 0;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    if (m_gas[i] == "CO2") {
      withCO2 = true;
      cCO2 = m_fraction[i];
      iCO2 = i;
    } else if (m_gas[i] == "CH4") {
      withCH4 = true;
      cCH4 = m_fraction[i];
      iCH4 = i;
    } else if (m_gas[i] == "C2H6") {
      withC2H6 = true;
      cC2H6 = m_fraction[i];
      iC2H6 = i;
    } else if (m_gas[i] == "C2H2") {
      withC2H2 = true;
      cC2H2 = m_fraction[i];
      iC2H2 = i;
    } else if (m_gas[i] == "CF4") {
      withCF4 = true;
      cCF4 = m_fraction[i];
      iCF4 = i;
    } else if (m_gas[i] == "iC4H10") {
      withIso = true;
      cIso = m_fraction[i];
      iIso = i;
    }
  }

  if (withAr && withCO2) {
    // Partial density of CO2
    const double nQ = GetNumberDensity() * cCO2;
    for (int j = nDeexcitations; j--;) {
      std::string level = m_deexcitations[j].label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      if (!optData.GetPhotoabsorptionCrossSection(
               "CO2", m_deexcitations[j].energy, pacs, eta)) {
        pacs = eta = 0.;
      }
      const double pPenningWK = pow(eta, 2. / 5.);
      if (level == "Ar_1S5") {
        // Rate constant from Velazco et al., J. Chem. Phys. 69 (1978)
        const double kQ = 5.3e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_1S4") {
        // Rate constant from Velazco et al., J. Chem. Phys. 69 (1978)
        const double kQ = 5.0e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_1S3") {
        const double kQ = 5.9e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_1S2") {
        const double kQ = 7.4e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2P8") {
        // Rate constant from Sadeghi et al., J. Chem. Phys. 115 (2001)
        const double kQ = 6.4e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2P6") {
        // Rate constant from Sadeghi et al.
        const double kQ = 6.1e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2P5") {
        // Rate constant from Sadeghi et al.
        const double kQ = 6.6e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2P1") {
        // Rate constant from Sadeghi et al.
        const double kQ = 6.2e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        // Average of 4p rate constants from Sadeghi et al.
        const double kQ = 6.33e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (m_deexcitations[j].osc > 0.) {
        // Higher resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double m1 = ElectronMassGramme / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMassGramme / (m_rgas[iCO2] - 1.);
        // Compute the reduced mass.
        double mR = m1 * m2 / (m1 + m2);
        mR /= AtomicMassUnit;
        const double uA =
            (RydbergEnergy / m_deexcitations[j].energy) * m_deexcitations[j].osc;
        const double uQ =
            (2 * RydbergEnergy / m_deexcitations[j].energy) * pacs /
            (4 * Pi2 * FineStructureConstant * BohrRadius * BohrRadius);
        const double kQ =
            2.591e-19 * pow(uA * uQ, 2. / 5.) * pow(m_temperature / mR, 3. / 10.);
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by CO2 (W-K formula):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        double pPenning = pPenningWK;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        // Collision radii
        const double rAr3d = 436.e-10;
        const double rCO2 = 165.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr3d + rCO2, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iCO2] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = fit3dQCO2 * sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by CO2 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        double pPenning = fit3dEtaCO2;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        // Collision radii
        const double rAr5s = 635.e-10;
        const double rCO2 = 165.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr5s + rCO2, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iCO2] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = fit3dQCO2 * sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by CO2 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        double pPenning = fit3dEtaCO2;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      }
    }
  }
  if (withAr && withCH4) {
    // Partial density of methane
    const double nQ = GetNumberDensity() * cCH4;
    for (int j = nDeexcitations; j--;) {
      std::string level = m_deexcitations[j].label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      if (!optData.GetPhotoabsorptionCrossSection(
               "CH4", m_deexcitations[j].energy, pacs, eta)) {
        pacs = eta = 0.;
      }
      const double pPenningWK = pow(eta, 2. / 5.);
      if (level == "Ar_1S5") {
        // Rate constant from Chen and Setser, J. Phys. Chem. 95 (1991)
        const double kQ = 4.55e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_1S4") {
        // Rate constant from Velazco et al., J. Chem. Phys. 69 (1978)
        const double kQ = 4.5e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_1S3") {
        // Rate constant from Chen and Setser
        const double kQ = 5.30e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_1S2") {
        // Rate constant from Velazco et al.
        const double kQ = 5.7e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2P8") {
        // Rate constant from Sadeghi et al., J. Chem. Phys. 115 (2001)
        const double kQ = 7.4e-19;
        double pPenning = pPenningWK;
        if (pPenning > 0.) pPenning = fit4pEtaCH4;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P6") {
        // Rate constant from Sadeghi et al.
        const double kQ = 3.4e-19;
        double pPenning = pPenningWK;
        if (pPenning > 0.) pPenning = fit4pEtaCH4;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P5") {
        // Rate constant from Sadeghi et al.
        const double kQ = 6.0e-19;
        double pPenning = pPenningWK;
        if (pPenning > 0.) pPenning = fit4pEtaCH4;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P1") {
        // Rate constant from Sadeghi et al.
        const double kQ = 9.3e-19;
        double pPenning = pPenningWK;
        if (pPenning > 0.) pPenning = fit4pEtaCH4;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        // Average of rate constants given by Sadeghi et al.
        const double kQ = 6.53e-19;
        double pPenning = pPenningWK;
        if (pPenning > 0.) pPenning = fit4pEtaCH4;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (m_deexcitations[j].osc > 0.) {
        // Higher resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double m1 = ElectronMassGramme / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMassGramme / (m_rgas[iCH4] - 1.);
        // Compute the reduced mass.
        double mR = m1 * m2 / (m1 + m2);
        mR /= AtomicMassUnit;
        const double uA =
            (RydbergEnergy / m_deexcitations[j].energy) * m_deexcitations[j].osc;
        const double uQ =
            (2 * RydbergEnergy / m_deexcitations[j].energy) * pacs /
            (4 * Pi2 * FineStructureConstant * BohrRadius * BohrRadius);
        const double kQ =
            2.591e-19 * pow(uA * uQ, 2. / 5.) * pow(m_temperature / mR, 3. / 10.);
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by CH4 (W-K formula):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        double pPenning = pPenningWK;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        // Collision radii
        const double rAr3d = 436.e-10;
        const double rCH4 = 190.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr3d + rCH4, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iCH4] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = fit3dQCH4 * sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by CH4 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        double pPenning = fit3dEtaCH4;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        // Collision radii
        const double rAr5s = 635.e-10;
        const double rCH4 = 190.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr5s + rCH4, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iCH4] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = fit3dQCH4 * sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by CH4 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        double pPenning = fit3dEtaCH4;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      }
    }
  }
  if (withAr && withC2H6) {
    // Partial density of ethane
    const double nQ = GetNumberDensity() * cC2H6;
    for (int j = nDeexcitations; j--;) {
      std::string level = m_deexcitations[j].label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      if (!optData.GetPhotoabsorptionCrossSection(
               "C2H6", m_deexcitations[j].energy, pacs, eta)) {
        pacs = eta = 0.;
      }
      const double pPenningWK = pow(eta, 2. / 5.);
      if (level == "Ar_1S5") {
        // Rate constant from Chen and Setser, J. Phys. Chem. 95 (1991)
        const double kQ = 5.29e-19;
        const double pPenning = pPenningWK;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_1S4") {
        // Rate constant from Velazco et al., J. Chem. Phys. 69 (1978)
        const double kQ = 6.2e-19;
        const double pPenning = pPenningWK;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_1S3") {
        // Rate constant from Chen and Setser
        const double kQ = 6.53e-19;
        const double pPenning = fit4sEtaC2H6;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_1S2") {
        // Rate constant from Velazco et al.
        const double kQ = 10.7e-19;
        const double pPenning = pPenningWK;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P8") {
        // Rate constant from Sadeghi et al., J. Chem. Phys. 115 (2001)
        const double kQ = 9.2e-19;
        double pPenning = fit4pEtaC2H6;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P6") {
        // Rate constant from Sadeghi et al.
        const double kQ = 4.8e-19;
        double pPenning = fit4pEtaC2H6;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P5") {
        // Rate constant from Sadeghi et al.
        const double kQ = 9.9e-19;
        double pPenning = fit4pEtaC2H6;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P1") {
        // Rate constant from Sadeghi et al.
        const double kQ = 11.0e-19;
        double pPenning = fit4pEtaC2H6;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        // Average of rate constants given by Sadeghi et al.
        const double kQ = 8.7e-19;
        double pPenning = fit4pEtaC2H6;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (m_deexcitations[j].osc > 0.) {
        // Higher resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double m1 = ElectronMassGramme / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMassGramme / (m_rgas[iC2H6] - 1.);
        // Compute the reduced mass.
        double mR = m1 * m2 / (m1 + m2);
        mR /= AtomicMassUnit;
        const double uA =
            (RydbergEnergy / m_deexcitations[j].energy) * m_deexcitations[j].osc;
        const double uQ =
            (2 * RydbergEnergy / m_deexcitations[j].energy) * pacs /
            (4 * Pi2 * FineStructureConstant * BohrRadius * BohrRadius);
        const double kQ =
            2.591e-19 * pow(uA * uQ, 2. / 5.) * pow(m_temperature / mR, 3. / 10.);
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by C2H6 (W-K formula):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        // Collision radii
        const double rAr3d = 436.e-10;
        const double rC2H6 = 195.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr3d + rC2H6, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iC2H6] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = fit3dQC2H6 * sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by C2H6 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        double pPenning = fit3dEtaC2H6;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        // Collision radii
        const double rAr5s = 635.e-10;
        const double rC2H6 = 195.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr5s + rC2H6, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iC2H6] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = fit3dQC2H6 * sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by C2H6 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        double pPenning = fit3dEtaC2H6;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      }
    }
  }
  if (withAr && withIso) {
    // Partial density of isobutane
    const double nQ = GetNumberDensity() * cIso;
    for (int j = nDeexcitations; j--;) {
      std::string level = m_deexcitations[j].label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      // Use n-butane as approximation for isobutane.
      if (!optData.GetPhotoabsorptionCrossSection(
               "nC4H10", m_deexcitations[j].energy, pacs, eta)) {
        pacs = eta = 0.;
      }
      const double pPenningWK = pow(eta, 2. / 5.);
      if (level == "Ar_1S5") {
        // Rate constant from
        // Piper et al., J. Chem. Phys. 59 (1973), 3323-3340
        const double kQ = 7.1e-19;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_1S4") {
        // Rate constant from Piper et al.
        const double kQ = 6.1e-19;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_1S3") {
        // Rate constant for n-butane from
        // Velazco et al., J. Chem. Phys. 69 (1978)
        const double kQ = 8.5e-19;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_1S2") {
        // Rate constant from Piper et al.
        const double kQ = 11.0e-19;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P8") {
        // Rate constant for ethane
        const double kEth = 9.2e-19;
        // Ar radius [pm]
        const double r4p = 340.;
        // Molecular radii [pm]
        const double rEth = 195.;
        const double rIso = 250.;
        // Masses [amu]
        const double mAr = 39.9;
        const double mEth = 30.1;
        const double mIso = 58.1;
        // Estimate rate constant for isobutane.
        double kQ = kEth;
        kQ *= pow((r4p + rIso) / (r4p + rEth), 2) *
              sqrt((mEth / mIso) * (mAr + mIso) / (mAr + mEth));
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Estim. rate constant for coll. deexcitation of\n"
                    << "    " << level << " by iC4H10:\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P6") {
        // Rate constant for ethane
        const double kEth = 4.8e-19;
        // Ar radius [pm]
        const double r4p = 340.;
        // Molecular radii [pm]
        const double rEth = 195.;
        const double rIso = 250.;
        // Masses [amu]
        const double mAr = 39.9;
        const double mEth = 30.1;
        const double mIso = 58.1;
        // Estimate rate constant for isobutane.
        double kQ = kEth;
        kQ *= pow((r4p + rIso) / (r4p + rEth), 2) *
              sqrt((mEth / mIso) * (mAr + mIso) / (mAr + mEth));
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Estim. rate constant for coll. deexcitation of\n"
                    << "    " << level << " by iC4H10:\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P5") {
        // Rate constant for ethane
        const double kEth = 9.9e-19;
        // Ar radius [pm]
        const double r4p = 340.;
        // Molecular radii [pm]
        const double rEth = 195.;
        const double rIso = 250.;
        // Masses [amu]
        const double mAr = 39.9;
        const double mEth = 30.1;
        const double mIso = 58.1;
        // Estimate rate constant for isobutane.
        double kQ = kEth;
        kQ *= pow((r4p + rIso) / (r4p + rEth), 2) *
              sqrt((mEth / mIso) * (mAr + mIso) / (mAr + mEth));
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Estim. rate constant for coll. deexcitation of\n"
                    << "    " << level << " by iC4H10:\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P1") {
        // Rate constant for Ethane
        const double kEth = 11.0e-19;
        // Ar radius [pm]
        const double r4p = 340.;
        // Molecular radii [pm]
        const double rEth = 195.;
        const double rIso = 250.;
        // Masses [amu]
        const double mAr = 39.9;
        const double mEth = 30.1;
        const double mIso = 58.1;
        // Estimate rate constant for isobutane.
        double kQ = kEth;
        kQ *= pow((r4p + rIso) / (r4p + rEth), 2) *
              sqrt((mEth / mIso) * (mAr + mIso) / (mAr + mEth));
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Estim. rate constant for coll. deexcitation of\n"
                    << "    " << level << " by iC4H10:\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        // Rate constante for ethane
        const double kEth = 5.5e-19;
        // Ar radius [pm]
        const double r4p = 340.;
        // Molecular radii [pm]
        const double rEth = 195.;
        const double rIso = 250.;
        // Masses [amu]
        const double mAr = 39.9;
        const double mEth = 30.1;
        const double mIso = 58.1;
        // Estimate rate constant for isobutane.
        double kQ = kEth;
        kQ *= pow((r4p + rIso) / (r4p + rEth), 2) *
              sqrt((mEth / mIso) * (mAr + mIso) / (mAr + mEth));
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Estim. rate constant for coll. deexcitation of\n"
                    << "    " << level << " by iC4H10:\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (m_deexcitations[j].osc > 0.) {
        // Higher resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double m1 = ElectronMassGramme / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMassGramme / (m_rgas[iIso] - 1.);
        // Compute the reduced mass.
        double mR = m1 * m2 / (m1 + m2);
        mR /= AtomicMassUnit;
        const double uA =
            (RydbergEnergy / m_deexcitations[j].energy) * m_deexcitations[j].osc;
        const double uQ =
            (2 * RydbergEnergy / m_deexcitations[j].energy) * pacs /
            (4 * Pi2 * FineStructureConstant * BohrRadius * BohrRadius);
        const double kQ =
            2.591e-19 * pow(uA * uQ, 2. / 5.) * pow(m_temperature / mR, 3. / 10.);
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by C4H10 (W-K formula):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        // Collision radii
        const double rAr3d = 436.e-10;
        const double rIso = 250.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr3d + rIso, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iIso] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by iC4H10 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        // Collision radii
        const double rAr5s = 635.e-10;
        const double rIso = 250.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr5s + rIso, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iIso] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by iC4H10 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      }
    }
  }
  if (withAr && withC2H2) {
    // Partial density of acetylene
    const double nQ = GetNumberDensity() * cC2H2;
    for (int j = nDeexcitations; j--;) {
      std::string level = m_deexcitations[j].label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      if (!optData.GetPhotoabsorptionCrossSection(
               "C2H2", m_deexcitations[j].energy, pacs, eta)) {
        pacs = eta = 0.;
      }
      const double pPenningWK = pow(eta, 2. / 5.);
      if (level == "Ar_1S5") {
        // Rate constant from Velazco et al., J. Chem. Phys. 69 (1978)
        const double kQ = 5.6e-19;
        // Branching ratio for ionization according to
        // Jones et al., J. Phys. Chem. 89 (1985)
        // p = 0.61, p = 0.74 (agrees roughly with WK estimate)
        const double pPenning = 0.61;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_1S4") {
        // Rate constant from Velazco et al.
        const double kQ = 4.6e-19;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_1S3") {
        const double kQ = 5.6e-19;
        const double pPenning = 0.61;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_1S2") {
        // Rate constant from Velazco et al.
        const double kQ = 8.7e-19;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P8") {
        // Rate constant from Sadeghi et al., J. Chem. Phys. 115 (2001)
        const double kQ = 5.0e-19;
        const double pPenning = 0.3;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P6") {
        // Rate constant from Sadeghi et al.
        const double kQ = 5.7e-19;
        const double pPenning = 0.3;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P5") {
        // Rate constant from Sadeghi et al.
        const double kQ = 6.0e-19;
        const double pPenning = 0.3;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P1") {
        // Rate constant from Sadeghi et al.
        const double kQ = 5.3e-19;
        const double pPenning = 0.3;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        // Average of rate constants given by Sadeghi et al.
        const double kQ = 5.5e-19;
        const double pPenning = 0.3;
        m_deexcitations[j].p.push_back(kQ * nQ * pPenning);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenning));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (m_deexcitations[j].osc > 0.) {
        // Higher resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double m1 = ElectronMassGramme / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMassGramme / (m_rgas[iC2H2] - 1.);
        // Compute the reduced mass.
        double mR = m1 * m2 / (m1 + m2);
        mR /= AtomicMassUnit;
        const double uA =
            (RydbergEnergy / m_deexcitations[j].energy) * m_deexcitations[j].osc;
        const double uQ =
            (2 * RydbergEnergy / m_deexcitations[j].energy) * pacs /
            (4 * Pi2 * FineStructureConstant * BohrRadius * BohrRadius);
        const double kQ =
            2.591e-19 * pow(uA * uQ, 2. / 5.) * pow(m_temperature / mR, 3. / 10.);
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by C2H2 (W-K formula):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        // Collision radii
        const double rAr3d = 436.e-10;
        const double rC2H2 = 165.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr3d + rC2H2, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iC2H2] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by C2H2 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        // Collision radii
        const double rAr5s = 635.e-10;
        const double rC2H2 = 165.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr5s + rC2H2, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iC2H2] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by C2H2 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ * pPenningWK);
        m_deexcitations[j].p.push_back(kQ * nQ * (1. - pPenningWK));
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollIon);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 2;
      }
    }
  }
  if (withAr && withCF4) {
    // Partial density of CF4
    const double nQ = GetNumberDensity() * cCF4;
    for (int j = nDeexcitations; j--;) {
      std::string level = m_deexcitations[j].label;
      // Photoabsorption cross-section and ionisation yield
      double pacs = 0., eta = 0.;
      if (!optData.GetPhotoabsorptionCrossSection(
               "CF4", m_deexcitations[j].energy, pacs, eta)) {
        pacs = eta = 0.;
      }
      if (level == "Ar_1S5") {
        // Rate constant from Chen and Setser
        const double kQ = 0.33e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_1S3") {
        // Rate constant from Chen and Setser
        const double kQ = 0.26e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2P8") {
        // Rate constant from Sadeghi et al.
        const double kQ = 1.7e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2P6") {
        // Rate constant from Sadeghi et al.
        const double kQ = 1.7e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2P5") {
        // Rate constant from Sadeghi et al.
        const double kQ = 1.6e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2P1") {
        // Rate constant from Sadeghi et al.
        const double kQ = 2.2e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2P10" || level == "Ar_2P9" || level == "Ar_2P7" ||
                 level == "Ar_2P4" || level == "Ar_2P3" || level == "Ar_2P2") {
        // Average of 4p rate constants from Sadeghi et al.
        const double kQ = 1.8e-19;
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (m_deexcitations[j].osc > 0.) {
        // Resonance levels
        // Calculate rate constant from Watanabe-Katsuura formula.
        const double m1 = ElectronMassGramme / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMassGramme / (m_rgas[iCF4] - 1.);
        // Compute the reduced mass.
        double mR = m1 * m2 / (m1 + m2);
        mR /= AtomicMassUnit;
        const double uA =
            (RydbergEnergy / m_deexcitations[j].energy) * m_deexcitations[j].osc;
        const double uQ =
            (2 * RydbergEnergy / m_deexcitations[j].energy) * pacs /
            (4 * Pi2 * FineStructureConstant * BohrRadius * BohrRadius);
        const double kQ =
            2.591e-19 * pow(uA * uQ, 2. / 5.) * pow(m_temperature / mR, 3. / 10.);
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by CF4 (W-K formula):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_3D6" || level == "Ar_3D3" || level == "Ar_3D4!" ||
                 level == "Ar_3D4" || level == "Ar_3D1!!" ||
                 level == "Ar_3D1!" || level == "Ar_3S1!!!!" ||
                 level == "Ar_3S1!!" || level == "Ar_3S1!!!") {
        // Non-resonant 3d levels
        // Collision radii
        const double rAr3d = 436.e-10;
        const double rCF4 = 235.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr3d + rCF4, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iCF4] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by CF4 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      } else if (level == "Ar_2S5" || level == "Ar_2S3") {
        // Non-resonant 5s levels
        // Collision radii
        const double rAr5s = 635.e-10;
        const double rCF4 = 190.e-10;
        // Hard sphere cross-section
        const double sigma = pow(rAr5s + rCF4, 2) * Pi;
        // Reduced mass
        const double m1 = ElectronMass / (m_rgas[iAr] - 1.);
        const double m2 = ElectronMass / (m_rgas[iCF4] - 1.);
        const double mR = m1 * m2 / (m1 + m2);
        // Relative velocity
        const double vel = SpeedOfLight * sqrt(8. * BoltzmannConstant *
                                               m_temperature / (Pi * mR));
        const double kQ = sigma * vel;
        if (m_debug) {
          std::cout << m_className << "::ComputeDeexcitationTable:\n";
          std::cout << "    Rate constant for coll. deexcitation of\n"
                    << "    " << level << " by CF4 (hard sphere):\n"
                    << "      " << kQ << " cm3 ns-1\n";
        }
        m_deexcitations[j].p.push_back(kQ * nQ);
        m_deexcitations[j].final.push_back(-1);
        m_deexcitations[j].type.push_back(DxcTypeCollNonIon);
        m_deexcitations[j].nChannels += 1;
      }
    }
  }

  if ((m_debug || verbose) && nDeexcitations > 0) {
    std::cout << m_className << "::ComputeDeexcitationTable:\n";
    std::cout << "      Level  Energy [eV]   "
              << "                 Lifetimes [ns]\n";
    std::cout << "                          "
              << "  Total    Radiative       "
              << "     Collisional\n";
    std::cout << "                               "
              << "                Ionisation  Transfer      Loss\n";
  }

  for (unsigned int i = 0; i < nDeexcitations; ++i) {
    // Calculate the total decay rate of each level.
    m_deexcitations[i].rate = 0.;
    double fRad = 0.;
    double fCollIon = 0., fCollTransfer = 0., fCollLoss = 0.;
    for (int j = m_deexcitations[i].nChannels; j--;) {
      m_deexcitations[i].rate += m_deexcitations[i].p[j];
      if (m_deexcitations[i].type[j] == DxcTypeRad) {
        fRad += m_deexcitations[i].p[j];
      } else if (m_deexcitations[i].type[j] == DxcTypeCollIon) {
        fCollIon += m_deexcitations[i].p[j];
      } else if (m_deexcitations[i].type[j] == DxcTypeCollNonIon) {
        if (m_deexcitations[i].final[j] < 0) {
          fCollLoss += m_deexcitations[i].p[j];
        } else {
          fCollTransfer += m_deexcitations[i].p[j];
        }
      } else {
        std::cerr << m_className << "::ComputeDeexcitationTable:\n";
        std::cerr << "    Unknown type of deexcitation channel (level "
                  << m_deexcitations[i].label << ")\n";
        std::cerr << "    Program bug!\n";
      }
    }
    if (m_deexcitations[i].rate > 0.) {
      // Print the radiative and collisional decay rates.
      if (m_debug || verbose) {
        std::cout << std::setw(12) << m_deexcitations[i].label << "  "
                  << std::fixed << std::setprecision(3) << std::setw(7)
                  << m_deexcitations[i].energy << "  " << std::setw(10)
                  << 1. / m_deexcitations[i].rate << "  ";
        if (fRad > 0.) {
          std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                    << 1. / fRad << " ";
        } else {
          std::cout << "---------- ";
        }
        if (fCollIon > 0.) {
          std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                    << 1. / fCollIon << " ";
        } else {
          std::cout << "---------- ";
        }
        if (fCollTransfer > 0.) {
          std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                    << 1. / fCollTransfer << " ";
        } else {
          std::cout << "---------- ";
        }
        if (fCollLoss > 0.) {
          std::cout << std::fixed << std::setprecision(3) << std::setw(10)
                    << 1. / fCollLoss << "\n";
        } else {
          std::cout << "---------- \n";
        }
      }
      // Normalise the decay branching ratios.
      for (int j = 0; j < m_deexcitations[i].nChannels; ++j) {
        m_deexcitations[i].p[j] /= m_deexcitations[i].rate;
        if (j > 0) m_deexcitations[i].p[j] += m_deexcitations[i].p[j - 1];
      }
    }
  }
}

void MediumMagboltz::ComputeDeexcitation(int iLevel, int& fLevel) {

  if (!m_useDeexcitation) {
    std::cerr << m_className << "::ComputeDeexcitation:\n";
    std::cerr << "    Deexcitation is disabled.\n";
    return;
  }

  // Make sure that the tables are updated.
  if (m_isChanged) {
    if (!Mixer()) {
      std::cerr << m_className << "::ComputeDeexcitation:\n";
      std::cerr << "    Error calculating the collision rates table.\n";
      return;
    }
    m_isChanged = false;
  }

  if (iLevel < 0 || iLevel >= (int)m_nTerms) {
    std::cerr << m_className << "::ComputeDeexcitation:\n";
    std::cerr << "    Level index is out of range.\n";
    return;
  }

  iLevel = m_iDeexcitation[iLevel];
  if (iLevel < 0 || iLevel >= (int)m_deexcitations.size()) {
    std::cerr << m_className << "::ComputeDeexcitation:\n";
    std::cerr << "    Level is not deexcitable.\n";
    return;
  }

  ComputeDeexcitationInternal(iLevel, fLevel);
  if (fLevel >= 0 && fLevel < (int)m_deexcitations.size()) {
    fLevel = m_deexcitations[fLevel].level;
  }
}

void MediumMagboltz::ComputeDeexcitationInternal(int iLevel, int& fLevel) {

  nDeexcitationProducts = 0;
  m_dxcProducts.clear();

  dxcProd newDxcProd;
  newDxcProd.s = 0.;
  newDxcProd.t = 0.;

  fLevel = iLevel;
  const int nDeexcitations = m_deexcitations.size();
  while (iLevel >= 0 && iLevel < nDeexcitations) {
    if (m_deexcitations[iLevel].rate <= 0. ||
        m_deexcitations[iLevel].nChannels <= 0) {
      // This level is a dead end.
      fLevel = iLevel;
      return;
    }
    // Determine the de-excitation time.
    newDxcProd.t += -log(RndmUniformPos()) / m_deexcitations[iLevel].rate;
    // Select the transition.
    fLevel = -1;
    int type = DxcTypeRad;
    const double r = RndmUniform();
    for (int j = 0; j < m_deexcitations[iLevel].nChannels; ++j) {
      if (r <= m_deexcitations[iLevel].p[j]) {
        fLevel = m_deexcitations[iLevel].final[j];
        type = m_deexcitations[iLevel].type[j];
        break;
      }
    }
    if (type == DxcTypeRad) {
      // Radiative decay
      newDxcProd.type = DxcProdTypePhoton;
      newDxcProd.energy = m_deexcitations[iLevel].energy;
      if (fLevel >= 0) {
        // Decay to a lower lying excited state.
        newDxcProd.energy -= m_deexcitations[fLevel].energy;
        if (newDxcProd.energy < Small) newDxcProd.energy = Small;
        m_dxcProducts.push_back(newDxcProd);
        ++nDeexcitationProducts;
        // Proceed with the next level in the cascade.
        iLevel = fLevel;
      } else {
        // Decay to ground state.
        double delta = RndmVoigt(0., m_deexcitations[iLevel].sDoppler,
                                 m_deexcitations[iLevel].gPressure);
        while (newDxcProd.energy + delta < Small ||
               fabs(delta) >= m_deexcitations[iLevel].width) {
          delta = RndmVoigt(0., m_deexcitations[iLevel].sDoppler,
                            m_deexcitations[iLevel].gPressure);
        }
        newDxcProd.energy += delta;
        m_dxcProducts.push_back(newDxcProd);
        ++nDeexcitationProducts;
        // Deexcitation cascade is over.
        fLevel = iLevel;
        return;
      }
    } else if (type == DxcTypeCollIon) {
      // Ionisation electron
      newDxcProd.type = DxcProdTypeElectron;
      newDxcProd.energy = m_deexcitations[iLevel].energy;
      if (fLevel >= 0) {
        // Associative ionisation
        newDxcProd.energy -= m_deexcitations[fLevel].energy;
        if (newDxcProd.energy < Small) newDxcProd.energy = Small;
        ++m_nPenning;
        m_dxcProducts.push_back(newDxcProd);
        ++nDeexcitationProducts;
        // Proceed with the next level in the cascade.
        iLevel = fLevel;
      } else {
        // Penning ionisation
        newDxcProd.energy -= m_minIonPot;
        if (newDxcProd.energy < Small) newDxcProd.energy = Small;
        ++m_nPenning;
        m_dxcProducts.push_back(newDxcProd);
        ++nDeexcitationProducts;
        // Deexcitation cascade is over.
        fLevel = iLevel;
        return;
      }
    } else if (type == DxcTypeCollNonIon) {
      // Proceed with the next level in the cascade.
      iLevel = fLevel;
    } else {
      std::cerr << m_className << "::ComputeDeexcitationInternal:\n";
      std::cerr << "    Unknown deexcitation channel type (" << type << ").\n";
      std::cerr << "    Program bug!\n";
      // Abort the deexcitation calculation.
      fLevel = iLevel;
      return;
    }
  }
}

bool MediumMagboltz::ComputePhotonCollisionTable(const bool verbose) {

  OpticalData data;
  double cs;
  double eta;

  // Atomic density
  const double dens = GetNumberDensity();

  // Reset the collision rate arrays.
  m_cfTotGamma.clear();
  m_cfTotGamma.resize(nEnergyStepsGamma, 0.);
  m_cfGamma.clear();
  m_cfGamma.resize(nEnergyStepsGamma);
  for (int j = nEnergyStepsGamma; j--;) m_cfGamma[j].clear();
  csTypeGamma.clear();

  nPhotonTerms = 0;
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    const double prefactor = dens * SpeedOfLight * m_fraction[i];
    // Check if optical data for this gas is available.
    std::string gasname = m_gas[i];
    if (gasname == "iC4H10") {
      gasname = "nC4H10";
      if (m_debug || verbose) {
        std::cout << m_className << "::ComputePhotonCollisionTable:\n";
        std::cout << "    Photoabsorption cross-section for "
                  << "iC4H10 not available.\n";
        std::cout << "    Using n-butane cross-section instead.\n";
      }
    }
    if (!data.IsAvailable(gasname)) return false;
    csTypeGamma.push_back(i * nCsTypesGamma + PhotonCollisionTypeIonisation);
    csTypeGamma.push_back(i * nCsTypesGamma + PhotonCollisionTypeInelastic);
    nPhotonTerms += 2;
    for (int j = 0; j < nEnergyStepsGamma; ++j) {
      // Retrieve total photoabsorption cross-section and ionisation yield.
      data.GetPhotoabsorptionCrossSection(gasname, (j + 0.5) * m_eStepGamma, cs,
                                          eta);
      m_cfTotGamma[j] += cs * prefactor;
      // Ionisation
      m_cfGamma[j].push_back(cs * prefactor * eta);
      // Inelastic absorption
      m_cfGamma[j].push_back(cs * prefactor * (1. - eta));
    }
  }

  // If requested, write the cross-sections to file.
  if (m_useCsOutput) {
    std::ofstream csfile;
    csfile.open("csgamma.txt", std::ios::out);
    for (int j = 0; j < nEnergyStepsGamma; ++j) {
      csfile << (j + 0.5) * m_eStepGamma << "  ";
      for (int i = 0; i < nPhotonTerms; ++i) csfile << m_cfGamma[j][i] << "  ";
      csfile << "\n";
    }
    csfile.close();
  }

  // Calculate the cumulative rates.
  for (int j = 0; j < nEnergyStepsGamma; ++j) {
    for (int i = 0; i < nPhotonTerms; ++i) {
      if (i > 0) m_cfGamma[j][i] += m_cfGamma[j][i - 1];
    }
  }

  if (m_debug || verbose) {
    std::cout << m_className << "::ComputePhotonCollisionTable:\n";
    std::cout << "    Energy [eV]      Mean free path [um]\n";
    for (int i = 0; i < 10; ++i) {
      const double imfp =
          m_cfTotGamma[(2 * i + 1) * nEnergyStepsGamma / 20] / SpeedOfLight;
      std::cout << "    " << std::fixed << std::setw(10) << std::setprecision(2)
                << (2 * i + 1) * m_eFinalGamma / 20 << "    " << std::setw(18)
                << std::setprecision(4);
      if (imfp > 0.) {
        std::cout << 1.e4 / imfp << "\n";
      } else {
        std::cout << "------------\n";
      }
    }
    std::cout << std::resetiosflags(std::ios_base::floatfield);
  }

  if (!m_useDeexcitation) return true;

  // Conversion factor from oscillator strength to cross-section
  const double f2cs =
      FineStructureConstant * 2 * Pi2 * HbarC * HbarC / ElectronMass;
  // Discrete absorption lines
  int nResonanceLines = 0;
  const unsigned int nDeexcitations = m_deexcitations.size();
  for (unsigned int i = 0; i < nDeexcitations; ++i) {
    if (m_deexcitations[i].osc < Small) continue;
    const double prefactor =
        dens * SpeedOfLight * m_fraction[m_deexcitations[i].gas];
    m_deexcitations[i].cf = prefactor * f2cs * m_deexcitations[i].osc;
    // Compute the line width due to Doppler broadening.
    const double mgas = ElectronMass / (m_rgas[m_deexcitations[i].gas] - 1.);
    const double wDoppler = sqrt(BoltzmannConstant * m_temperature / mgas);
    m_deexcitations[i].sDoppler = wDoppler * m_deexcitations[i].energy;
    // Compute the half width at half maximum due to resonance broadening.
    //   A. W. Ali and H. R. Griem, Phys. Rev. 140, 1044
    //   A. W. Ali and H. R. Griem, Phys. Rev. 144, 366
    const double kResBroad = 1.92 * Pi * sqrt(1. / 3.);
    m_deexcitations[i].gPressure = kResBroad * FineStructureConstant *
                                 pow(HbarC, 3) * m_deexcitations[i].osc * dens *
                                 m_fraction[m_deexcitations[i].gas] /
                                 (ElectronMass * m_deexcitations[i].energy);
    // Make an estimate for the width within which a photon can be
    // absorbed by the line
    // const int nWidths = 1000;
    const double nWidths = fitLineCut;
    // Calculate the FWHM of the Voigt distribution according to the
    // approximation formula given in
    // Olivero and Longbothum, J. Quant. Spectr. Rad. Trans. 17, 233-236
    const double fwhmGauss = m_deexcitations[i].sDoppler * sqrt(2. * log(2.));
    const double fwhmLorentz = m_deexcitations[i].gPressure;
    const double fwhmVoigt =
        0.5 * (1.0692 * fwhmLorentz + sqrt(0.86639 * fwhmLorentz * fwhmLorentz +
                                           4 * fwhmGauss * fwhmGauss));
    m_deexcitations[i].width = nWidths * fwhmVoigt;
    ++nResonanceLines;
  }

  if (nResonanceLines <= 0) {
    std::cerr << m_className << "::ComputePhotonCollisionTable:\n";
    std::cerr << "    No resonance lines found.\n";
    return true;
  }

  if (m_debug || verbose) {
    std::cout << m_className << "::ComputePhotonCollisionTable:\n";
    std::cout << "    Discrete absorption lines:\n";
    std::cout << "      Energy [eV]        Line width (FWHM) [eV]  "
              << "    Mean free path [um]\n";
    std::cout << "                            Doppler    Pressure   "
              << "   (peak)     \n";
    for (unsigned int i = 0; i < nDeexcitations; ++i) {
      if (m_deexcitations[i].osc < Small) continue;
      const double imfpP = (m_deexcitations[i].cf / SpeedOfLight) *
                           TMath::Voigt(0., m_deexcitations[i].sDoppler,
                                        2 * m_deexcitations[i].gPressure);
      std::cout << "      " << std::fixed << std::setw(6)
                << std::setprecision(3) << m_deexcitations[i].energy << " +/- "
                << std::scientific << std::setprecision(1)
                << m_deexcitations[i].width << "   " << std::setprecision(2)
                << 2 * sqrt(2 * log(2.)) * m_deexcitations[i].sDoppler << "   "
                << std::scientific << std::setprecision(3)
                << 2 * m_deexcitations[i].gPressure << "  " << std::fixed
                << std::setw(10) << std::setprecision(4);
      if (imfpP > 0.) {
        std::cout << 1.e4 / imfpP;
      } else {
        std::cout << "----------";
      }
      std::cout << "\n";
    }
  }

  return true;
}

void MediumMagboltz::RunMagboltz(const double e, const double bmag,
                                 const double btheta, const int ncoll,
                                 bool verbose, double& vx, double& vy,
                                 double& vz, double& dl, double& dt,
                                 double& alpha, double& eta, double& lor,
                                 double& vxerr, double& vyerr, double& vzerr, 
                                 double& dlerr, double& dterr, 
                                 double& alphaerr, double& etaerr, 
                                 double& lorerr, double& alphatof) {

  // Initialize the values.
  vx = vy = vz = 0.;
  dl = dt = 0.;
  alpha = eta = alphatof = 0.;
  lor = 0.;
  vxerr = vyerr = vzerr = 0.;
  dlerr = dterr = 0.;
  alphaerr = etaerr = 0.;
  lorerr = 0.;

  // Set input parameters in Magboltz common blocks.
  Magboltz::inpt_.nGas = m_nComponents;
  Magboltz::inpt_.nStep = 4000;
  Magboltz::inpt_.nAniso = 2;

  Magboltz::inpt_.tempc = m_temperature - ZeroCelsius;
  Magboltz::inpt_.torr = m_pressure;
  Magboltz::inpt_.ipen = 0;
  Magboltz::setp_.nmax = ncoll;

  Magboltz::setp_.efield = e;
  // Convert from Tesla to kGauss.
  Magboltz::bfld_.bmag = bmag * 10.;
  // Convert from radians to degree.
  Magboltz::bfld_.btheta = btheta * 180. / Pi;

  // Set the gas composition in Magboltz.
  for (unsigned int i = 0; i < m_nComponents; ++i) {
    int ng = 0;
    if (!GetGasNumberMagboltz(m_gas[i], ng)) {
      std::cerr << m_className << "::RunMagboltz:\n"
                << "    Gas " << m_gas[i] << " has no corresponding"
                << " gas number in Magboltz.\n";
      return;
    }
    Magboltz::gasn_.ngasn[i] = ng;
    Magboltz::ratio_.frac[i] = 100 * m_fraction[i];
  }

  // Call Magboltz internal setup routine.
  Magboltz::setup1_();

  // Calculate the max. energy in the table.
  if (e * m_temperature / (293.15 * m_pressure) > 15) {
    // If E/p > 15 start with 8 eV.
    Magboltz::inpt_.efinal = 8.;
  } else {
    Magboltz::inpt_.efinal = 0.5;
  }
  Magboltz::setp_.estart = Magboltz::inpt_.efinal / 50.;

  long long ielow = 1;
  while (ielow == 1) {
    Magboltz::mixer_();
    if (bmag == 0. || btheta == 0. || fabs(btheta) == Pi) {
      Magboltz::elimit_(&ielow);
    } else if (btheta == HalfPi) {
      Magboltz::elimitb_(&ielow);
    } else {
      Magboltz::elimitc_(&ielow);
    }
    if (ielow == 1) {
      // Increase the max. energy.
      Magboltz::inpt_.efinal *= sqrt(2.);
      Magboltz::setp_.estart = Magboltz::inpt_.efinal / 50.;
    }
  }

  if (m_debug || verbose) Magboltz::prnter_();

  // Run the Monte Carlo calculation.
  if (bmag == 0.) {
    Magboltz::monte_();
  } else if (btheta == 0. || btheta == Pi) {
    Magboltz::montea_();
  } else if (btheta == HalfPi) {
    Magboltz::monteb_();
  } else {
    Magboltz::montec_();
  }
  if (m_debug || verbose) Magboltz::output_();

  // If attachment or ionisation rate is greater than sstmin,
  // include spatial gradients in the solution.
  const double sstmin = 30.;
  const double epscale = 760. * m_temperature / (m_pressure * 293.15);
  double alpp = Magboltz::ctowns_.alpha * epscale;
  double attp = Magboltz::ctowns_.att * epscale;
  bool useSST = false;
  if (fabs(alpp - attp) > sstmin || alpp > sstmin || attp > sstmin) {
    useSST = true;
    if (bmag == 0.) {
      Magboltz::alpcalc_();
    } else if (btheta == 0. || btheta == Pi) {
      Magboltz::alpclca_();
    } else if (btheta == HalfPi) {
      Magboltz::alpclcb_();
    } else {
      Magboltz::alpclcc_();
    }
    // Calculate the (effective) TOF Townsend coefficient.
    double alphapt = Magboltz::tofout_.ralpha;
    double etapt = Magboltz::tofout_.rattof;
    double fc1 =
        1.e5 * Magboltz::tofout_.tofwr / (2. * Magboltz::tofout_.tofdl);
    double fc2 = 1.e12 * (alphapt - etapt) / Magboltz::tofout_.tofdl;
    alphatof = fc1 - sqrt(fc1 * fc1 - fc2);
  }
  if (m_debug || verbose) Magboltz::output2_();

  // Velocities. Convert to cm / ns.
  vx = Magboltz::vel_.wx * 1.e-9;
  vxerr = Magboltz::velerr_.dwx;
  vy = Magboltz::vel_.wy * 1.e-9;
  vyerr = Magboltz::velerr_.dwy;
  vz = Magboltz::vel_.wz * 1.e-9;
  vzerr = Magboltz::velerr_.dwz;

  // Calculate the Lorentz angle.
  const double forcalc = vx * vx + vy * vy;
  double elvel = sqrt(forcalc + vz * vz);
  if (forcalc != 0 && elvel != 0) {
    lor = acos(vz / elvel);
    const double ainlorerr = sqrt(forcalc * forcalc * vzerr * vzerr + 
                                  vx * vx * vx * vx * vxerr * vxerr + 
                                  vy * vy * vy * vy * vyerr * vyerr);
    lorerr = vz * ainlorerr/ elvel / elvel / sqrt (forcalc) / lor;
  }

  // Diffusion coefficients.
  // dt = sqrt(0.2 * Magboltz::difvel_.diftr / vz) * 1.e-4;
  dt = sqrt(0.2 * 0.5 * (Magboltz::diflab_.difxx + Magboltz::diflab_.difyy) / 
            vz) * 1.e-4;
  dterr = Magboltz::diferl_.dfter;
  // dl = sqrt(0.2 * Magboltz::difvel_.difln / vz) * 1.e-4;
  dl = sqrt(0.2 * Magboltz::diflab_.difzz / vz) * 1.e-4;
  dlerr = Magboltz::diferl_.dfler;
  // Diffusion tensor.
  // SBOL(1)=2D-6*DIFZZ*TORR/VBOL
  // SBOL(2)=2D-6*DIFXX*TORR/VBOL
  // SBOL(3)=2D-6*DIFYY*TORR/VBOL
  // SBOL(4)=2D-6*DIFXZ*TORR/VBOL
  // SBOL(5)=2D-6*DIFYZ*TORR/VBOL
  // SBOL(6)=2D-6*DIFXY*TORR/VBOL
  alpha = Magboltz::ctowns_.alpha;
  alphaerr = Magboltz::ctwner_.alper;
  eta = Magboltz::ctowns_.att;
  etaerr = Magboltz::ctwner_.atter;

  // Print the results.
  if (m_debug) {
    std::cout << m_className << "::RunMagboltz:\n    Results:\n";
    std::cout << "      Drift velocity along E:           " << std::right
              << std::setw(10) << std::setprecision(6) << vz << " cm/ns +/- "
              << std::setprecision(2) << vzerr << "%\n";
    std::cout << "      Drift velocity along Bt:          " << std::right
              << std::setw(10) << std::setprecision(6) << vx << " cm/ns +/- "
              << std::setprecision(2) << vxerr << "%\n";
    std::cout << "      Drift velocity along ExB:         " << std::right
              << std::setw(10) << std::setprecision(6) << vy << " cm/ns +/- "
              << std::setprecision(2) << vyerr << "%\n";
    std::cout << "      Longitudinal diffusion:           " << std::right
              << std::setw(10) << std::setprecision(6) << dl << " cm1/2 +/- "
              << std::setprecision(2) << dlerr << "%\n";
    std::cout << "      Transverse diffusion:             " << std::right
              << std::setw(10) << std::setprecision(6) << dt << " cm1/2 +/- "
              << std::setprecision(2) << dterr << "%\n";
    std::cout << "      Lorentz Angle:           " << std::right
              << std::setw(10) << std::setprecision(6) << (lor / Pi * 180.) 
              << " degree  +/- " << std::setprecision(2) << lorerr << "%\n";
    if (useSST) {
      std::cout << "      Townsend coefficient (SST):       " << std::right
                << std::setw(10) << std::setprecision(6) << alpha
                << " cm-1  +/- " << std::setprecision(2) << alphaerr << "%\n";
      std::cout << "      Attachment coefficient (SST):     " << std::right
                << std::setw(10) << std::setprecision(6) << eta << " cm-1  +/- "
                << std::setprecision(2) << etaerr << "%\n";
      std::cout << "      Eff. Townsend coefficient (TOF):  " << std::right
                << std::setw(10) << std::setprecision(6) << alphatof
                << " cm-1\n";
    } else {
      std::cout << "      Townsend coefficient:             " << std::right
                << std::setw(10) << std::setprecision(6) << alpha
                << " cm-1  +/- " << std::setprecision(2) << alphaerr << "%\n";
      std::cout << "      Attachment coefficient:           " << std::right
                << std::setw(10) << std::setprecision(6) << eta << " cm-1  +/- "
                << std::setprecision(2) << etaerr << "%\n";
    }
  }
}

void MediumMagboltz::GenerateGasTable(const int numColl, const bool verbose) {

  // Set the reference pressure and temperature.
  m_pressureTable = m_pressure;
  m_temperatureTable = m_temperature;

  // Initialize the parameter arrays.
  const unsigned int nEfields = m_eFields.size();
  const unsigned int nBfields = m_bFields.size();
  const unsigned int nAngles = m_bAngles.size();
  InitParamArrays(nEfields, nBfields, nAngles, tabElectronVelocityE, 0.);
  InitParamArrays(nEfields, nBfields, nAngles, tabElectronVelocityB, 0.);
  InitParamArrays(nEfields, nBfields, nAngles, tabElectronVelocityExB, 0.);
  InitParamArrays(nEfields, nBfields, nAngles, tabElectronDiffLong, 0.);
  InitParamArrays(nEfields, nBfields, nAngles, tabElectronDiffTrans, 0.);
  InitParamArrays(nEfields, nBfields, nAngles, tabElectronLorentzAngle, 0.);
  InitParamArrays(nEfields, nBfields, nAngles, tabElectronTownsend, -30.);
  InitParamArrays(nEfields, nBfields, nAngles, m_tabTownsendNoPenning, -30.);
  InitParamArrays(nEfields, nBfields, nAngles, tabElectronAttachment, -30.);

  m_hasElectronVelocityE = true;
  m_hasElectronVelocityB = true;
  m_hasElectronVelocityExB = true;
  m_hasElectronDiffLong = true;
  m_hasElectronDiffTrans = true;
  m_hasElectronAttachment = true;
  m_hasElectronLorentzAngle = true;

  m_hasExcRates = false;
  m_tabExcRates.clear();
  m_excitationList.clear();
  m_hasIonRates = false;
  m_tabIonRates.clear();
  m_ionisationList.clear();

  m_hasIonMobility = false;
  m_hasIonDissociation = false;
  m_hasIonDiffLong = false;
  m_hasIonDiffTrans = false;

  // gasBits = "TFTTFTFTTTFFFFFF";
  // The version number is 11 because there are slight
  // differences between the way these gas files are written
  // and the ones from Garfield. This is mainly in the way
  // the gas tables are stored.
  // versionNumber = 11;

  double vx = 0., vy = 0., vz = 0.;
  double difl = 0., dift = 0.;
  double alpha = 0., eta = 0.;
  double lor = 0.;
  double vxerr = 0., vyerr = 0., vzerr = 0.;
  double diflerr = 0., difterr = 0.;
  double alphaerr = 0., etaerr = 0.;
  double alphatof = 0.;
  double lorerr = 0.;

  // Run through the grid of E- and B-fields and angles.
  for (unsigned int i = 0; i < nEfields; ++i) {
    for (unsigned int j = 0; j < nAngles; ++j) {
      for (unsigned int k = 0; k < nBfields; ++k) {
        if (m_debug) {
          std::cout << m_className << "::GenerateGasTable:\n"
                    << "    E = " << m_eFields[i] << " V/cm, B = " 
                    << m_bFields[k] << " T, angle: " << m_bAngles[j] << " rad\n";
        }
        RunMagboltz(m_eFields[i], m_bFields[k], m_bAngles[j], numColl, verbose, vx,
                    vy, vz, difl, dift, alpha, eta, lor, vxerr, vyerr, vzerr,
                    diflerr, difterr, alphaerr, etaerr, lorerr, alphatof);
        tabElectronVelocityE[j][k][i] = vz;
        tabElectronVelocityExB[j][k][i] = vy;
        tabElectronVelocityB[j][k][i] = vx;
        tabElectronDiffLong[j][k][i] = difl;
        tabElectronDiffTrans[j][k][i] = dift;
        tabElectronLorentzAngle[j][k][i] = lor;
        if (alpha > 0.) {
          tabElectronTownsend[j][k][i] = log(alpha);
          m_tabTownsendNoPenning[j][k][i] = log(alpha);
        } else {
          tabElectronTownsend[j][k][i] = -30.;
          m_tabTownsendNoPenning[j][k][i] = -30.;
        }
        if (eta > 0.) {
          tabElectronAttachment[j][k][i] = log(eta);
        } else {
          tabElectronAttachment[j][k][i] = -30.;
        }
      }
    }
  }
}
}
