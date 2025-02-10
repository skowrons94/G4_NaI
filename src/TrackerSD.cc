#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4EventManager.hh"
#include "G4AnalysisManager.hh"
#include <vector>
#include <queue>

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

TrackerSD::TrackerSD(const G4String& name, const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
    collectionName.insert(hitsCollectionName);
}

TrackerSD::~TrackerSD() 
{}

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 
    G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(hcID, fHitsCollection);
    
    // Clear steps at the start of each event
    fSteps.clear();
}

G4bool TrackerSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{  
    G4double edep = step->GetTotalEnergyDeposit();
    
    // Create hit object
    TrackerHit* newHit = new TrackerHit();
    newHit->SetTrackID(step->GetTrack()->GetTrackID());
    newHit->SetChamberNb(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());
    newHit->SetEdep(edep);
    newHit->SetPos(step->GetPostStepPoint()->GetPosition());
    fHitsCollection->insert(newHit);
    
    // Store step information
    StepInfo info;
    info.trackID = step->GetTrack()->GetTrackID();
    info.parentID = step->GetTrack()->GetParentID();
    info.edep = edep;
    info.particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    
    fSteps.push_back(info);
    
    return true;
}

G4double TrackerSD::CalculateTotalEnergy(G4int currentID, 
                                       const std::map<G4int, std::vector<G4int>>& daughters,
                                       const std::map<G4int, G4double>& energyDeposits) {
    G4double totalEnergy = 0.0;
    
    // Add energy deposit of current particle
    auto energyIt = energyDeposits.find(currentID);
    if (energyIt != energyDeposits.end()) {
        totalEnergy += energyIt->second;
    }
    
    // Recursively add energy of all daughters
    auto daughterIt = daughters.find(currentID);
    if (daughterIt != daughters.end()) {
        for (G4int daughterID : daughterIt->second) {
            totalEnergy += CalculateTotalEnergy(daughterID, daughters, energyDeposits);
        }
    }
    
    return totalEnergy;
}

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();
    
    // Get max track ID
    G4int maxTrackID = 0;
    for (const auto& step : fSteps) {
        maxTrackID = std::max(maxTrackID, step.trackID);
    }

    // Create a lookup map for steps by trackID for quick access
    std::unordered_map<G4int, const StepInfo*> trackMap;
    for (const auto& step : fSteps) {
        trackMap[step.trackID] = &step;
    }

    // Set of analyzed track IDs to avoid re-analysis
    std::unordered_set<G4int> analyzedTracks;

    // Recursive function to accumulate edep for a given trackID's secondaries
    std::function<G4double(G4int)> accumulateEdep = [&](G4int trackID) -> G4double {
        G4double totalEdep = 0.0;
        
       // Find all the dughters of the current track
       // And do it recursively until all daughters are found
        std::queue<G4int> toVisit;
        toVisit.push(trackID);
        while (!toVisit.empty()) {
            G4int currentID = toVisit.front();
            toVisit.pop();
            
            // Find the daughetrs
            // So when the current ID is the parent ID of another track
            // Then we add the energy deposition of the daughter to the total energy deposition
            for (const auto& step : fSteps) {
                // Check if in already analyzed
                if (analyzedTracks.find(step.trackID) != analyzedTracks.end()) continue;
                // Else check if the current ID is the parent ID of the step
                if (step.parentID == currentID) {
                    toVisit.push(step.trackID);
                    analyzedTracks.insert(step.trackID);
                    totalEdep += step.edep;
                }
            }
        }
        return totalEdep;
    };

    // Main loop to go through each track ID
    for (G4int trackID = 0; trackID <= maxTrackID; ++trackID) {
        // Check if the track is already analyzed
        if (analyzedTracks.find(trackID) != analyzedTracks.end()) continue;

        // Check if the trackID exists in trackMap
        if (trackMap.find(trackID) != trackMap.end()) {
            const auto& step = trackMap[trackID];
            
            // Check if particle is gamma or alpha
            if (step->particleName == "gamma" || step->particleName == "alpha") {
                analyzedTracks.insert(trackID);  // Mark the particle as analyzed

                // Accumulate energy deposition including all secondaries
                G4double totalEdep = step->edep + accumulateEdep(trackID);

                // Store the result in the corresponding histogram
                if (step->particleName == "gamma") {
                    analysis->FillH1( 12, totalEdep );
                } else if (step->particleName == "alpha") {
                    analysis->FillH1( 13, totalEdep );
                }
            }
        }
    }
    
}