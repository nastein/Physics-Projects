#include "g4PSIDetectorBase.hh"

g4PSIDetectorBase::g4PSIDetectorBase(G4String label) {
    label_ = label;
    //  log_ = NULL;
    collID_ = -1;
    mother_volume_ = NULL;
    detector_log_ = NULL;
    HitRequired_ = false;
    info_wiki_ = false;
}


void g4PSIDetectorBase::Info(std::ofstream *file, int n) {
    info_wiki_n_ = n;
    info_wiki_ = true;
    info_wiki_file_ = file;
    Info();
}

void g4PSIDetectorBase::InfoTitle(G4String name) {
    if (info_wiki_) {
        *info_wiki_file_ << "|| (" << info_wiki_n_ << ") " << name << " (" << label_ << ") || ||\n";
    } else {
        G4cout << "NAME: " << name << " (" << label_ << ")" << G4endl;
    }
}

void g4PSIDetectorBase::InfoParDouble(G4String par, G4double val, G4String unit) {
    if (info_wiki_) {
        *info_wiki_file_ << "| " << par << " | " << val << " " << unit << " |\n";
    } else {
        G4cout << par << " = " << val << " " << unit << G4endl;
    }
}

void g4PSIDetectorBase::InfoPar2Double(G4String par, G4double val1, G4double val2, G4String unit) {
    if (info_wiki_) {
        *info_wiki_file_ << "| " << par << " | " <<
        val1 << " " << unit << " x " <<
        val2 << " " << unit << " |\n";
    } else {
        G4cout << par << " = " <<
        val1 << " " << unit << " x " <<
        val2 << " " << unit <<
        G4endl;
    }
}

void g4PSIDetectorBase::InfoPar3Double(G4String par, G4double val1, G4double val2, G4double val3, G4String unit) {
    if (info_wiki_) {
        *info_wiki_file_ << "| " << par << " | " <<
        val1 << " " << unit << " x " <<
        val2 << " " << unit << " x " <<
        val3 << " " << unit << " |\n";
    } else {
        G4cout << par << " = " <<
        val1 << " " << unit << " x " <<
        val2 << " " << unit << " x " <<
        val3 << " " << unit <<
        G4endl;
    }
}

void g4PSIDetectorBase::InfoParBool(G4String par, G4bool val) {
    if (info_wiki_) {
        *info_wiki_file_ << "| " << par << " | " << (val ? "yes" : "no") << " |\n";
    } else {
        G4cout << par << " = " << (val ? "yes" : "no") << G4endl;
    }
}

void g4PSIDetectorBase::InfoParString(G4String par, G4String val) {
    if (info_wiki_) {
        *info_wiki_file_ << "| " << par << " | " << val << " |\n";
    } else {
        G4cout << par << " = " << val << G4endl;
    }
}

void g4PSIDetectorBase::InfoParInt(G4String par, G4int val, G4String unit) {
    if (info_wiki_) {
        *info_wiki_file_ << "| " << par << " | " << val << " " << unit << " |\n";
    } else {
        G4cout << par << " = " << val << " " << unit << G4endl;
    }
}

void g4PSIDetectorBase::InfoEnd() {
    if (info_wiki_) {
        *info_wiki_file_ << G4endl;
    }
}