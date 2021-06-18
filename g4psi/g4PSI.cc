/*! \mainpage Geant4 Simulation of the PSI muon-proton experiment
 *
 *
 * "This product includes software developed by Members of the Geant4
 *  Collaboration ( http://cern.ch/geant4 )."
 */


//#undef G4VIS_USE


#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4StepLimiterPhysics.hh"

#include "g4PSIDetectorConstruction.hh"
#include "g4PSIPrimaryGeneratorAction.hh"
#include "g4PSISteppingAction.hh"
#include "g4PSIRunAction.hh"
#include "g4PSIEventAction.hh"
#include "g4PSIActionInitialization.hh"
#include "g4PSIAnalysisManager.hh"

#include "FTFP_BERT_MUSE.hh"
#include "G4GenericBiasingPhysics.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

/// main for g4PSI Geant4 simulation code

// http://stackoverflow.com/questions/865668/parse-command-line-arguments
char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

// ---------------------------------------------------------------------------------------

int main(int argc,char** argv)
{
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    
    // Construct the default run manager
    G4RunManager* runManager = new G4RunManager;
    
    // set mandatory initialization classes
    runManager->SetUserInitialization(new g4PSIDetectorConstruction);
    
    FTFP_BERT_MUSE* physicsList = new FTFP_BERT_MUSE;
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());

    int argshift = 0;
    if(cmdOptionExists(argv, argv+argc, "--bias"))
    {
        argshift++;
        // augment physicsList with biasing facilities:
        G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
        // -- Create list of physics processes to be biased: only CoulombScat. in this case:
        std::vector<G4String> processToBias;
        processToBias.push_back("CoulombScat");
        //      processToBias.push_back("eBrem");
        //      processToBias.push_back("eIoni");

        // -- Pass the list to the G4GenericBiasingPhysics, which will wrap the CoulombScat
        // -- process of e- and e+ to activate the biasing of it:
        biasingPhysics->Bias("e-", processToBias);
        biasingPhysics->Bias("e+", processToBias);
        biasingPhysics->Bias("mu-", processToBias);
        biasingPhysics->Bias("mu+", processToBias);
        physicsList->RegisterPhysics(biasingPhysics);
        g4PSIAnalysisManager::getInstance()->SetBiasingPhysics(true, 1000.0);

        G4cout << "      ********************************************************* " << G4endl;
        G4cout << "      ********** processes are wrapped for biasing ************ " << G4endl;
        G4cout << "      ********************************************************* " << G4endl;
    }
    else
    {
        g4PSIAnalysisManager::getInstance()->SetBiasingPhysics(false, 1.0);
        G4cout << "      ************************************************* " << G4endl;
        G4cout << "      ********** processes are not wrapped ************ " << G4endl;
        G4cout << "      ************************************************* " << G4endl;
    }

    runManager->SetUserInitialization(physicsList);

    runManager->SetUserAction(new g4PSIPrimaryGeneratorAction());
    runManager->SetUserAction(new g4PSIRunAction());
    runManager->SetUserAction(new g4PSIEventAction());
    runManager->SetUserAction(new g4PSISteppingAction());


    // Initialize G4 kernel
    
    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    if(argc == 1 + argshift)  // Define (G)UI terminal for interactive mode
    {
#ifdef G4VIS_USE
        // Visualization, if you choose to have it!
        G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();
#endif
        
#ifdef G4UI_USE
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
        std::cout << "COMMAND: Settings.mac\n";
        UImanager->ApplyCommand("/control/execute Settings.mac");
#endif
        if (ui->IsGUI())
            std::cout << "COMMAND: gui.mac\n";
        UImanager->ApplyCommand("/control/execute gui.mac");
        ui->SessionStart();
        delete ui;
#endif
        
#ifdef G4VIS_USE
        G4cout << "ending visManager" << G4endl;
        delete visManager;
#endif
        
    }
    else   // Batch mode
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1 + argshift];
        UImanager->ApplyCommand(command+fileName);
    }
    
    G4cout << "ending runManager" << G4endl;
    delete runManager;
    G4cout.flush();
    return 0;
}
