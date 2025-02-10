// TrackerSD.hh
#ifndef TrackerSD_h
#define TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "HistoManager.hh"
#include "TrackerHit.hh"
#include <vector>
#include <map>
#include <set>

// Structure to store step information
struct StepInfo {
    G4int trackID;
    G4int parentID;
    G4double edep;
    G4String particleName;
};

class TrackerSD : public G4VSensitiveDetector
{
  public:
    TrackerSD(const G4String& name, const G4String& hitsCollectionName);
    virtual ~TrackerSD();
  
    virtual void Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    TrackerHitsCollection* fHitsCollection;
    HistoManager* fHistoManager;
    
    // Vector to store all steps in the event
    std::vector<StepInfo> fSteps;
    
    // Helper function to calculate total energy for a particle and its daughters
    G4double CalculateTotalEnergy(G4int currentID, const std::map<G4int, std::vector<G4int>>& daughters, 
                                 const std::map<G4int, G4double>& energyDeposits);
    
    // Helper function to find all daughter particles
    void FindDaughters(G4int currentID, std::set<G4int>& allDaughters, 
                      const std::map<G4int, std::vector<G4int>>& daughters);
};

#endif