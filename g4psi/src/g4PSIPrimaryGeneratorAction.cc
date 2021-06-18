///
///  \todo{http://hypernews.slac.stanford.edu/HyperNews/geant4/get/runmanage/264/1.html}
///

#include "g4PSIPrimaryGeneratorAction.hh"
#include "g4PSIPrimaryGeneratorMessenger.hh"
#include "g4PSIAnalysisManager.hh"
//#include "g4PSISCWall.hh"
//#include "g4PSIWC.hh"
#include "g4PSIDetectorParts.hh"

#include "G4Event.hh"
#include "G4KaonPlus.hh"
#include "G4PionPlus.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"

#include "TLorentzVector.h"

using namespace CLHEP;

g4PSIPrimaryGeneratorAction::g4PSIPrimaryGeneratorAction()
{
    gunMode_ = WALLUNIFORM;
    beamMixMode_ = DEFAULT_PARTICLE;
    max_nPrimaryParticle_ = 1;
    
    beamMomentumCenter_ = 153 * MeV;
    beamRelativeMomentumSpread_ = 0.00;
    path0_ = 24.5 * m; // Path length from channel opening to target center
    z0_beamline = -1.5*m; // Beamline mode only; Z origin of events
    proton_jitter = 10.0*ns; // Time width of proton pulse jitter
    target_radius_ = 3*cm; // Radius that the full4pi beam will be produced from
    
    event_seed1_ = 0;
    event_seed2_ = 0;
    load_seeds_ = false;
    close_seed_file_ = false;
    seed_filename_ = "seedfile.dat";
    turtle_filename_ = "";
    
    beam_box_x0_ = 0.0 * mm;
    beam_box_y0_ = 0.0 * mm;
    beam_box_dx_ = 12.0 * mm;
    beam_box_dy_ = 12.0 * mm;
    
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
    
    const G4int n_particle = 1;
    particleGun_ = new G4ParticleGun(n_particle);
    G4cout << "Test before " << particleGun_->GetParticleEnergy() << " " <<  particleGun_->GetParticleMomentum() << G4endl;
    particleGun_->SetParticleDefinition(particle);
    particleGun_->SetParticleEnergy(0);
    particleGun_->SetParticleMomentum(beamMomentumCenter_);
    particleGun_->SetParticlePosition(G4ThreeVector(0, 0, z0_beamline));
    particleGun_->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1.0));
    
    // create a messenger for this class
    gunMessenger = new g4PSIPrimaryGeneratorMessenger(this);
}


g4PSIPrimaryGeneratorAction::~g4PSIPrimaryGeneratorAction()
{
    delete particleGun_;
}


void g4PSIPrimaryGeneratorAction::GetBeamParticle(G4double &charge, G4double &mass, G4int &id) {
    
    G4ParticleDefinition* particle = NULL;
    
    if (beamMixMode_ == BEAM_MIX_PLUS ||
        beamMixMode_ == BEAM_MIX_MINUS) {
        
        G4String name;
        if (beamMixMode_ == BEAM_MIX_PLUS) {
            G4double choice = CLHEP::RandFlat::shoot(0.0, 3.0);
            if (choice < 1.00) {
                name = "e+";
            } else if (choice > 1.00 && choice < 2.00) {
                name = "mu+";
            } else {
                name = "pi+";
            };
        } else if (beamMixMode_ == BEAM_MIX_MINUS) {
            G4double choice = CLHEP::RandFlat::shoot(0.0, 3.0);
            if (choice < 1.00) {
                name = "e-";
            } else if (choice > 1.00 && choice < 2.00) {
                name = "mu-";
            } else {
                name = "pi-";
            }
        };
        
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        particle = particleTable->FindParticle(name);
        particleGun_->SetParticleDefinition(particle);
    } else {
        particle = particleGun_->GetParticleDefinition();
    };
    
    charge = particle->GetPDGCharge();
    mass =  particle->GetPDGMass();
    id = particle->GetPDGEncoding();
}


G4double g4PSIPrimaryGeneratorAction::GetBeamMomentum() {
    return beamMomentumCenter_ * (1.0 + CLHEP::RandFlat::shoot(-0.5,+0.5)*beamRelativeMomentumSpread_);
}


G4double g4PSIPrimaryGeneratorAction::GetRFTime(G4double mass, G4double pBeam, G4double z0) {
    const G4double c = 29.9792 * cm / ns;
    G4double E = sqrt(mass*mass + pBeam*pBeam);
    G4double beta = pBeam / E;
    G4double RF_time = (1/c) * ((path0_ + z0)/beta - path0_);
    
    // Optional: add jitter for proton pulse width
    if (proton_jitter > 0) {
        CLHEP::RandGauss::shoot();  // somehow we need to request a RandGauss in pairs, or our seed setting does not work.
        RF_time = RF_time + proton_jitter*CLHEP::RandGauss::shoot();
    }
    return RF_time;
}


void g4PSIPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    
    // clear the primary particle information.
    
    g4PSIAnalysisManager::getInstance()->ClearParticleGun();
    //std::cout << "Seeds before: " << "\n";
    //CLHEP::HepRandom::showEngineStatus();

    //long seed_holder[2];

    long seeds[2];
    seeds[0] = 0L;
    seeds[1] = 0L;
    if (load_seeds_) {
        if (close_seed_file_) {
            if (seed_file_.is_open()) seed_file_.close();
            close_seed_file_ = false;
        }
        if (!seed_file_.is_open()) {
            seed_file_.open(seed_filename_.c_str());
            if (seed_file_.is_open()) {
                std::cout << "USING SEED FILE: " << seed_filename_ << std::endl;
            } else {
                std::cout << "ERROR OPENING SEED FILE: " << seed_filename_ << std::endl;
            }
        }
        if (seed_file_.is_open()) {
            seed_file_ >> event_seed1_;
            seed_file_ >> event_seed2_;
            std::cout << "EVENT SEED: /gun/seeds " << event_seed1_ << " " << event_seed2_ << std::endl;
        }
        if (seed_file_.eof()) {
            std::cout << "ERROR: EOF\n";
        } else if (!seed_file_.is_open()) {
            std::cout << "ERROR: CAN NOT OPEN SEED FILE " << seed_filename_ << "\n";
        };
    };
    if (event_seed1_ > 0 || event_seed2_ > 0) {
        seeds[0] = event_seed1_;
        seeds[1] = event_seed2_;
        //G4cout << "before seeds[0] = " << seeds[0] << G4endl;
        //G4cout << "before Seeds[1] = " << seeds[1] << G4endl;
        //seed_holder[0] = seeds[0];
        //seed_holder[1] = seeds[1];
        CLHEP::HepRandom::setTheSeeds(seeds);
        //CLHEP::HepRandom::showEngineStatus();
        event_seed1_ = 0; // reset to zero
        event_seed2_ = 0; // reset to zero
    };
    //-----------------------------------------------------------------
    seeds[0] = *(CLHEP::HepRandom::getTheSeeds());
    seeds[1] = *(CLHEP::HepRandom::getTheSeeds()+1);
    G4int eventID = anEvent->GetEventID();
    g4PSIAnalysisManager::getInstance()->SetEventID(eventID, seeds[0], seeds[1]);
    
    for (int iPrimaryParticle = 0; iPrimaryParticle < max_nPrimaryParticle_; iPrimaryParticle++) {
        
        G4double mass = 0;
        G4double charge = 0;
        G4int particleID = 0;
        G4double x0 = 0;
        G4double y0 = 0;
        G4double z0 = 0;
        G4double RF_time = 0.0;
        G4double pin = GetBeamMomentum();
        G4double pout = pin;
        G4double beam_theta = 0.0;
        G4double beam_phi = 0.0;
        
        // Pick beam particle for the beam_mix modes
        // PDG IDs are: e+/e-    11+/-
        //              mu+/mu-  13+/-
        GetBeamParticle(charge, mass, particleID);
        
        if (gunMode_ == WALLUNIFORM ||
            gunMode_ == WALLUNIFORM_MOMDIST ||
            gunMode_ == FULL4PI ||
            gunMode_ == FULL4PI_MOMDIST ||
            gunMode_ == FULL4PI_FROM_CENTER ||
            gunMode_ == FULL4PI_DEFAULT ||
            gunMode_ == WALLUNIFORM_FROM_CENTER) {
            
            G4double px;
            G4double py;
            G4double pz;
            
            bool ok = false;
            while (!ok) {
                /// \todo for now the reaction vertex is in a rectangular volume, should be gaussian, possibly with divergence
                if (gunMode_ != WALLUNIFORM_FROM_CENTER &&
                    gunMode_ != FULL4PI_FROM_CENTER  &&
                    gunMode_ != FULL4PI_DEFAULT) {
                    bool in_target = false;
                    const double cut_h = 2.5 * cm;
                    const double cut_v = 3.5 * cm;
                    const double target_r = target_radius_;// was 3.0*cm
                    while (!(in_target && fabs(x0) < cut_h)) {
                        x0 = CLHEP::RandFlat::shoot(-target_r, target_r);
                        z0 = CLHEP::RandFlat::shoot(-target_r, target_r);
                        in_target = (x0*x0 + z0*z0) < target_r*target_r;
                    }
                    y0 = CLHEP::RandFlat::shoot(-cut_v,+cut_v);
                };
                // Pick spherically uniform momentum directions until the direction vector points into the wall
                //beam_theta = acos(CLHEP::RandFlat::shoot(-1,1));
                beam_theta = CLHEP::RandFlat::shoot(20., 100.);
                beam_phi = CLHEP::RandFlat::shoot(0.,360.);
                if (beam_phi > 180.0) beam_phi -= 360;
                beam_phi *= deg;
                beam_theta *= deg;
                px = std::sin(beam_theta)*std::cos(beam_phi);
                py = std::sin(beam_theta)*std::sin(beam_phi);
                pz = std::cos(beam_theta);
                ok = (gunMode_ == FULL4PI) || (gunMode_ == FULL4PI_MOMDIST) || (gunMode_ == FULL4PI_FROM_CENTER) || (gunMode_ == FULL4PI_DEFAULT) ||
                (beam_theta/deg > 20 && beam_theta/deg < 100 && fabs(beam_phi/deg) < 45.0);
            }
            
            if (gunMode_ != FULL4PI_DEFAULT) {
                particleGun_->SetParticlePosition(G4ThreeVector(x0, y0, z0));
            }
            particleGun_->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
            
            RF_time = GetRFTime(mass, pin, z0);
            
            if (gunMode_ == WALLUNIFORM_MOMDIST ||
                gunMode_ == FULL4PI_MOMDIST) {
                // shoot a flat momentum distribution with pin as maximum momentum.
                pout = CLHEP::RandFlat::shoot(0.0, pin);
            } else {
                pout = get_particle_momentum(pin, mass, beam_theta);
            };
            
        } else if (gunMode_ == TARGET_MUON_DECAY) {
            
            G4double px;
            G4double py;
            G4double pz;
            
            bool ok = false;
            
            while (!ok) {
                /// \todo for now the reaction vertex is in a rectangular volume, should be gaussian, possibly with divergence
                
                bool in_target = false;
                while (!in_target) {
                    const double rx = 2. * cm;
                    const double ry = 2. * cm;
                    x0 = CLHEP::RandFlat::shoot(-rx, rx);
                    y0 = CLHEP::RandFlat::shoot(-ry, ry);
                    in_target = ((x0/rx)*(x0/rx) + (y0/ry)*(y0/ry)) < 1.0;
                }
                z0 = CLHEP::RandFlat::shoot(-2.,+2.)*cm;
                
                // Pick spherically uniform momentum directions until the direction vector points into the wall
                beam_theta = acos(CLHEP::RandFlat::shoot(-1,1));
                beam_phi = CLHEP::RandFlat::shoot(0.,360.);
                if (beam_phi > 180.0) beam_phi -= 360;
                beam_phi *= deg;
                px = std::sin(beam_theta)*std::cos(beam_phi);
                py = std::sin(beam_theta)*std::sin(beam_phi);
                pz = std::cos(beam_theta);
                
                ok = (beam_theta > 20 * deg && beam_theta < 100 * deg && fabs(beam_phi) < 45.0 * deg);
            }
            
            particleGun_->SetParticlePosition(G4ThreeVector(x0, y0, z0));
            particleGun_->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
            
            RF_time = GetRFTime(mass, pin, z0);
            
            if (gunMode_ != WALLUNIFORM_MOMDIST) {
                pout = get_particle_momentum(pin, mass, beam_theta) * MeV;
            } else {
                // shoot a flat momentum distribution with pin as maximum momentum.
                pout = CLHEP::RandFlat::shoot(0.0, pin);
            };
            
        } else if (gunMode_ == COSMIC) {
            
            G4double px;
            G4double py;
            G4double pz;
            
            bool ok = false;
            while (!ok) {
                x0 = CLHEP::RandFlat::shoot(-5.,+5.)*m;
                z0 = CLHEP::RandFlat::shoot(-5.,+5.)*m;
                y0 = 3*m;
                beam_phi = CLHEP::RandFlat::shoot(0.,360. * deg);
                
                bool ok2 = false;
                while (!ok2) {
                    beam_theta = CLHEP::RandFlat::shoot(0.0, 90.0 * deg);
                    G4double theta_test = CLHEP::RandFlat::shoot(0.0, 1.0);
                    G4double ct = std::cos(beam_theta);
                    ok2 = theta_test < ct*ct;
                }
                px = std::sin(beam_theta)*std::cos(beam_phi);
                pz = std::sin(beam_theta)*std::sin(beam_phi);
                py = -std::cos(beam_theta);
                ok = true;
                //      ok = g4PSISCWall::getInstance()->hitsAnyWall(px, py, pz, x0, y0, z0);
            }
            
            particleGun_->SetParticlePosition(G4ThreeVector(x0, y0, z0));
            particleGun_->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
            
            pout = CLHEP::RandFlat::shoot(0.0, 2.0) * GeV;
            
        } else if (gunMode_ == WALLCENTER ||
                   gunMode_ == WALLCENTER_MOMDIST) {
            
            x0 = y0 = z0 = 0.0;
            G4double theta0 = 60*deg;  // debug, \todo, fix this later:  g4PSISCWall::getInstance()->GetSC_theta0();
            particleGun_->SetParticlePosition(G4ThreeVector(x0, y0, z0));
            G4ThreeVector v(sin(theta0), 0., cos(theta0));
            particleGun_->SetParticleMomentumDirection(v);
            if (gunMode_ == WALLCENTER_MOMDIST) {
                // shoot a flat momentum distribution with pin as maximum momentum.
                pout = CLHEP::RandFlat::shoot(0.0, pin);
            }
            
        } else if (gunMode_ == BEAMLINE_2015 ||
                   gunMode_ == BEAMLINE_2015_CENTER) {
            
            z0 = z0_beamline;
            
            G4double sig_xp = 0;
            G4double sig_yp = 0;
            G4double sig_x0 = 0;
            G4double sig_y0 = 0;
            G4double mean_xp = 0;
            G4double mean_yp = 0;
            G4double mean_x0 = 0;
            G4double mean_y0 = 0;
            G4double fpx =  0.;
            G4double fpy =  0.;
            
            if (charge <= 0) {
                
                // negatively charged particles, or neutrals, geantino
                
                if (pin < 130.0 * MeV) {
                    mean_xp = -8./1000.;
                    mean_yp = +8./1000.;
                    mean_x0 = 3. * mm;
                    mean_y0 = 2. * mm;
                    sig_xp = 18./1000.;
                    sig_yp =  6./1000.;
                    sig_x0 = 0. * mm;
                    sig_y0 = 6. * mm;
                    fpx =  0. * mm;
                    fpy =  -300. * mm;
                } else if (pin < 180.0 * MeV) {
                    mean_xp = -8./1000.;
                    mean_yp = +4./1000.;
                    mean_x0 = 2. * mm;
                    mean_y0 = 1. * mm;
                    sig_xp = 15./1000.;
                    sig_yp = 10./1000.;
                    sig_x0 = 2. * mm;
                    sig_y0 = 5. * mm;
                    fpx =  -20. * mm;
                    fpy =  -300. * mm;
                } else {
                    mean_xp = -2./1000.;
                    mean_yp = +4./1000.;
                    mean_x0 = 0. * mm;
                    mean_y0 = -1. * mm;
                    sig_xp = 10./1000.;
                    sig_yp = 7./1000.;
                    sig_x0 = 2. * mm;
                    sig_y0 = 4. * mm;
                    fpx =  -20. * mm;
                    fpy =  -400. * mm;
                }
            } else {
                
                // positively charged particles
                
                if (pin < 130.0 * MeV) {
                    mean_xp = -8./1000.;
                    mean_yp = +0./1000.;
                    mean_x0 = 1.5 * mm;
                    mean_y0 = 0. * mm;
                    sig_xp = 24./1000.;
                    sig_yp =  6./1000.;
                    sig_x0 = 2. * mm;
                    sig_y0 = 6. * mm;
                    fpx =  0. * mm;
                    fpy =  -300. * mm;
                } else if (pin < 180.0 * MeV) {
                    mean_xp = -8./1000.;
                    mean_yp = +0./1000.;
                    mean_x0 = 1. * mm;
                    mean_y0 = 0. * mm;
                    sig_xp = 16./1000.;
                    sig_yp = 10./1000.;
                    sig_x0 = 3. * mm;
                    sig_y0 = 5. * mm;
                    fpx =  -20. * mm;
                    fpy =  -300. * mm;
                } else {
                    mean_xp = +0./1000.;
                    mean_yp = +2./1000.;
                    mean_x0 = -1. * mm;
                    mean_y0 = -1. * mm;
                    sig_xp = 14./1000.;
                    sig_yp = 7./1000.;
                    sig_x0 = 2. * mm;
                    sig_y0 = 4. * mm;
                    fpx =  -20. * mm;
                    fpy =  -400. * mm;
                }
            }
            
            if (gunMode_ == BEAMLINE_2015_CENTER) {
                // no positional or directional offset
                mean_xp = 0.0;
                mean_yp = 0.0;
                mean_x0 = 0.0;
                mean_y0 = 0.0;
            }
            
            G4double xp = CLHEP::RandGauss::shoot(mean_xp, sig_xp);
            G4double yp = CLHEP::RandGauss::shoot(mean_yp, sig_yp);
            x0 = CLHEP::RandGauss::shoot(mean_x0, sig_x0) + (z0 - fpx) * xp;
            y0 = CLHEP::RandGauss::shoot(mean_y0, sig_y0) + (z0 - fpy) * yp;
            
            particleGun_->SetParticlePosition(G4ThreeVector(x0, y0, z0));
            particleGun_->SetParticleMomentumDirection(G4ThreeVector(xp, yp, 1.));
            
        } else if (gunMode_ == BEAMLINE_2014) {
            
            z0 = z0_beamline;
            
            G4double sig_xp = 0;
            G4double sig_yp = 0;
            G4double sig_x0 = 0;
            G4double sig_y0 = 0;
            G4double fpx =  0.;
            G4double fpy =  0.;
            
            if (pin < 130.0 * MeV) {
                sig_xp = 16./1000.;
                sig_yp = 10./1000.;
                sig_x0 = 0.2 * cm;
                sig_y0 = 0.7 * cm;
                fpx =  50. * cm;
                fpy =  30. * cm;
            } else if (pin < 180.0 * MeV) {
                sig_xp = 14./1000.;
                sig_yp = 10./1000.;
                sig_x0 = 0.2 * cm;
                sig_y0 = 0.8 * cm;
                fpx =  50. * cm;
                fpy =  0. * cm;
            } else {
                sig_xp = 12./1000.;
                sig_yp = 10./1000.;
                sig_x0 = 0.2 * cm;
                sig_y0 = 0.7 * cm;
                fpx =  30. * cm;
                fpy =  0. * cm;
            }
            
            G4double xp = CLHEP::RandGauss::shoot(0, sig_xp);
            G4double yp = CLHEP::RandGauss::shoot(0, sig_yp);
            x0 = CLHEP::RandGauss::shoot(0, sig_x0) + (z0 - fpx) * xp;
            y0 = CLHEP::RandGauss::shoot(0, sig_y0) + (z0 - fpy) * yp;
            
            particleGun_->SetParticlePosition(G4ThreeVector(x0, y0, z0));
            particleGun_->SetParticleMomentumDirection(G4ThreeVector(xp, yp, 1.));
            
        } else if (gunMode_ == BEAMLINE_PENCIL ||
                   gunMode_ == BEAMLINE_PENCIL_FLAT) {
            if (gunMode_ == BEAMLINE_PENCIL_FLAT) {
                pout = CLHEP::RandFlat::shoot(0.0, pin);
            }
            x0 = 0.;
            y0 = 0.;
            z0 = (particleGun_->GetParticlePosition()).z();
            particleGun_->SetParticlePosition(G4ThreeVector(x0, y0, z0));
            particleGun_->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
            
        } else if (gunMode_ == BEAMLINE_BOX) {
            x0 = beam_box_x0_ + 0.5 * CLHEP::RandFlat::shoot(-beam_box_dx_, beam_box_dx_);
            y0 = beam_box_y0_ + 0.5 * CLHEP::RandFlat::shoot(-beam_box_dy_, beam_box_dy_);
            z0 = (particleGun_->GetParticlePosition()).z();
            particleGun_->SetParticlePosition(G4ThreeVector(x0, y0, z0));
            particleGun_->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
        } else if (gunMode_ == TURTLE) {
            
            if (!turtle_file_.is_open()) {
                turtle_file_.open(turtle_filename_.c_str());
                if (turtle_file_.is_open()) {
                    std::cout << "USING TURTLE FILE: " << turtle_filename_ << std::endl;
                } else {
                    std::cout << "ERROR OPENING TURTLE FILE: " << turtle_filename_ << std::endl;
                }
            }
            if (turtle_file_.is_open()) {
                long n = 0;
                double w = 0;
                double xp = 0;
                double yp = 0;
                z0 = -140.0 * CLHEP::cm;
                turtle_file_ >> x0;
                turtle_file_ >> xp;
                turtle_file_ >> y0;
                turtle_file_ >> yp;
                turtle_file_ >> pout;
                turtle_file_ >> w;
                turtle_file_ >> n;
                
                x0 *= CLHEP::cm;
                y0 *= CLHEP::cm;
                xp *= CLHEP::mrad;
                yp *= CLHEP::mrad;
                pout *= CLHEP::MeV;
                
                double zfx = -x0 / xp + z0;
                double zfy = -y0 / yp + z0;
                
                std::cout << "TURTLE: " << x0 << " " << xp << " " << y0 << " " << yp << " " << pout << " f: " << zfx << " " << zfy << std::endl;

                particleGun_->SetParticlePosition(G4ThreeVector(x0, y0, z0));
                particleGun_->SetParticleMomentumDirection(G4ThreeVector(xp, yp, 1.));
            }
            if (turtle_file_.eof()) {
                std::cout << "ERROR: EOF\n";
            } else if (!turtle_file_.is_open()) {
                std::cout << "ERROR: CAN NOT OPEN TURTLE FILE " << turtle_filename_ << "\n";
            };
            
        }
    
        if (gunMode_ != DEFAULT_GUN) {
            //            RF_time = GetRFTime(mass, pin, z0);
            G4double kinEnergy = std::sqrt(pout*pout+mass*mass) - mass;
            particleGun_->SetParticleEnergy(kinEnergy);
        }
        
        g4PSIAnalysisManager::getInstance()->SetParticleGun(particleGun_);
        particleGun_->GeneratePrimaryVertex(anEvent);
    }
}

G4double g4PSIPrimaryGeneratorAction::get_particle_momentum(G4double pin, G4double mass_particle, G4double theta) {
    // theta in radian
    // from Ron:
    
    Double_t mass_p = 938.27;
    
    Double_t ein = TMath::Sqrt(TMath::Power(pin,2)+TMath::Power(mass_particle,2));
    
    Double_t xmssq  = TMath::Power(ein+mass_p,2) - TMath::Power(pin,2) + TMath::Power(mass_particle,2) - TMath::Power(mass_p,2);
    Double_t thr = theta;
    //calculate mu, p outgoing momenta and energies
    Double_t a = 4. * (TMath::Power(ein+mass_p,2) - TMath::Power(pin*TMath::Cos(thr),2));
    Double_t b = -4. * xmssq * pin * TMath::Cos(thr);
    Double_t c = 4. * TMath::Power(ein+mass_p,2) * TMath::Power(mass_particle,2) - TMath::Power(xmssq,2);
    Double_t bsqm4ac = TMath::Power(b,2) - 4.*a*c;
    Double_t pmuo = (-1.*b + TMath::Power(bsqm4ac,0.5)) / (2.*a);
    //  std::cout << "mass: " << mass_particle << " theta=" << theta/deg<< " p=" << pmuo << std::endl;
    return pmuo; // in MeV
}


void g4PSIPrimaryGeneratorAction::SetEventSeeds(G4long s1, G4long s2, bool ls, std::string fn) {
    event_seed1_ = s1; event_seed2_ = s2;
    load_seeds_ = ls;
    seed_filename_ = fn;
    close_seed_file_ = true;
}


void g4PSIPrimaryGeneratorAction::SetBeamMomentumCenter(G4double p) {
    beamMomentumCenter_ = p;
    particleGun_->SetParticleMomentum(p);
}


void g4PSIPrimaryGeneratorAction::SetBeamMomentumSpread(G4double dpp) {
    beamRelativeMomentumSpread_ = dpp;
}

void g4PSIPrimaryGeneratorAction::SetTargetRadius(G4double r) {
    target_radius_ = r;
}

void g4PSIPrimaryGeneratorAction::SetTurtleFile(G4String fn) {
    turtle_filename_ = fn;
}
