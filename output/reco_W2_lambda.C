#define reco_W2_lambda_cxx
#include "reco_W2_lambda.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <ROOT::Math::LorentzVector>

float reco_W2_lambda::GetInvariantMass(TVector3 p1, TVector3 p2) {
    // Assuming the tracks are pions and protons, reconstruct the invariant mass of the lambda
    float massPion = 0.13957; // pion mass in GeV/c^2
    float massProton = 0.93827; // proton mass in GeV/c^2

    TVector3 pionP3 = p1;
    TVector3 protonP3 = p2;
    double pionPMag = p1.Mag();
    double protonPMag = p2.Mag();
    float pionEnergy = std::sqrt(massPion*massPion + pionPMag*pionPMag);
    float protonEnergy = std::sqrt(massProton*massProton + protonPMag*protonPMag);
    LorentzVector protonP4;
    LorentzVector pionP4;

   pionP4.SetPxPyPzE(pionP3.X(), pionP3.Y(), pionP3.Z(), pionEnergy);
   protonP4.SetPxPyPzE(protonP3.X(), protonP3.Y(), protonP3.Z(), protonEnergy);
    
    LorentzVector lambdaP4 = protonP4 + pionP4;
    double lambdaInvariantMass = lambdaP4.M();

    return lambdaInvariantMass;
}

void reco_W2_lambda::Loop()
{
//   In a ROOT session, you can do:
//      root> .L reco_W2_lambda.C
//      root> reco_W2_lambda t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   // Define histograms

   TH1F* hTrueLambdaMassDirect = new TH1F("hTrueLambdaMassDirect", "True Lambda Invariant Mass (Direct); Mass (GeV/c^2); Entries", 100, 1.0, 1.2);
   TH1F* hTrueLambdaMassIndirect = new TH1F("hTrueLambdaMassIndirect", "True Lambda Invariant Mass (Decay Products); Mass (GeV/c^2); Entries", 100, 1.0, 1.2);
   TH1F* hRecoLambdaMass = new TH1F("hRecoLambdaMass", "Reconstructed Lambda Invariant Mass; Mass (GeV/c^2); Entries", 100, 1.0, 1.2);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // ------------------------------------------
      // True Lambda Invariant Mass
      // ------------------------------------------

      if(sampleType==0){ // only calculate true lambda invariant mass for signal events
   
      // method 1: directly from lambda
      for (size_t i = 0; i < truePDG->size(); ++i){
         if (truePDG->at(i) == 3122){ // Lambda PDG code
            TVector3 lambdaMomentum(truePx->at(i), truePy->at(i), truePz->at(i));
            double lambdaEnergy = trueE->at(i);
            double lambdaMassSquared = lambdaEnergy*lambdaEnergy - lambdaMomentum.Mag()*lambdaMomentum.Mag();
            double lambdaMassDirect = (lambdaMassSquared > 0) ? std::sqrt(lambdaMassSquared) : 0.0;
            hTrueLambdaMassDirect->Fill(lambdaMassDirect);
         }
      }

      // method 2: reconstruct from decay products (proton and pion)
      TVector3 protonMomentum(-9999., -9999., -9999.);
      TVector3 pionMomentum(-9999., -9999., -9999.);
      double protonEnergy;
      double pionEnergy;
      double protonMass;
      double pionMass;
      for (size_t i = 0; i < truePDG->size(); ++i){
         if (truePDG->at(i) == 2212){ // Proton PDG code
            protonMomentum.SetXYZ(truePx->at(i), truePy->at(i), truePz->at(i));
            protonEnergy = trueE->at(i);
            protonMass = trueMass->at(i);
         }
         else if (truePDG->at(i) == 211){ // Pion+ PDG code
            pionMomentum.SetXYZ(truePx->at(i), truePy->at(i), truePz->at(i));
            pionEnergy = trueE->at(i);
            pionMass = trueMass->at(i);
         }
      }

      if (protonMomentum.Mag() > 0 && pionMomentum.Mag() > 0){
         LorentzVector protonP4;
         LorentzVector pionP4;
         protonP4.SetPxPyPzE(protonMomentum.X(), protonMomentum.Y(), protonMomentum.Z(), protonEnergy);
         pionP4.SetPxPyPzE(pionMomentum.X(), pionMomentum.Y(), pionMomentum.Z(), pionEnergy);
         LorentzVector lambdaP4 = protonP4 + pionP4;
         double lambdaMassIndirect = lambdaP4.M();
         hTrueLambdaMassIndirect->Fill(lambdaMassIndirect);
      }

      } // end of sampleType==0


      // ------------------------------------------
      // Reco Invariant Mass
      // ------------------------------------------

      // Identify longest track in event -> designate as candidate muon
      int longestTrackID = -1;
      float longestTrackLength = -1.0;
      int longestTrackIdx = -1;
      for (size_t i = 0; i < pfpTrackScore->size(); ++i){
         if(pfpHasTrack && pfpTrackScore->at(i) > 0.5){
            if(pfpTrackLength->at(i) > longestTrackLength){
               longestTrackLength = pfpTrackLength->at(i);
               longestTrackID = pfpTrackID->at(i);
               longestTrackIdx = i;
            }
         }
      }

      int muonCandidateIdx = longestTrackIdx;

      // Designate other two tracks as proton and pion candidates
      int pionCandidateIdx = -1;
      int protonCandidateIdx = -1;
      float pionCandidateLength = -1.0;

      for (size_t i = 0; i < pfpTrackScore->size(); ++i){
         if(pfpHasTrack && pfpTrackScore->at(i) > 0.5 && i != muonCandidateIdx){
            if(pfpTrackLength->at(i) > pionCandidateLength){
               pionCandidateLength = pfpTrackLength->at(i);
               pionCandidateIdx = i;
            }

            else if(pfpTrackLength->at(i) <= pionCandidateLength && i != muonCandidateIdx){
               protonCandidateIdx = i;
            }
         }
      }

   } // End loop over TTree entries

   TCanvas* c1 = new TCanvas("c1", "Invariant Mass Histograms", 800, 600);
   c1->Divide(2, 1);
   c1->cd(1);
   hTrueLambdaMassDirect->SetLineColor(kBlue);
   hTrueLambdaMassDirect->Draw();
   c1->cd(2);
   hTrueLambdaMassIndirect->SetLineColor(kRed);
   hTrueLambdaMassIndirect->Draw();

   c1->SaveAs("lambda_invariant_mass_histograms.png");
}
