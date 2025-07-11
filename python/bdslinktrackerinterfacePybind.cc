//
// Created by Stewart Boogert on 07/06/2025.
//

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// #define CACHE_PTRS 1

#include "BDSLinkTrackerInterface.hh"
#include "BDSIMLink.hh"
#include "BDSLinkBunch.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "BDSException.hh"
#include "G4Electron.hh"
#include <BDSSamplerCustom.hh>

template <typename T>
T* make_ptr(py::array_t<T> &arr) {
  auto buf = arr.request();
  T* ptr = static_cast<T*>(buf.ptr);
  return ptr;
}

template <typename T>
inline void set_element(py::array_t<T> &arr, int index, T value) {
  T* ptr = make_ptr<T>(arr);
  ptr[index] = value;
}

void selectCollimator(BDSIMLink* link, const std::string& collimatorName)
{
    // This doesn't throw an error if the element doesn't exist
    std::cout << "Selecting collimator " << collimatorName << std::endl;
    link->SelectLinkElement(collimatorName);

    // Check if the element exists by querying the index: -1 means it doesn't exist
    if (link->GetLinkIndex(collimatorName) == -1)
        {throw std::runtime_error("Element not found " + collimatorName);}
}

BDSParticleDefinition* PrepareBDSParticleDefinition(long long int pdgIDIn, double momentumIn, 
                                                  double kineticEnergyIn, double ionChargeIn)
{
    G4int pdgID = (G4int) pdgIDIn;

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDefGeant;

    // Wrap in our class that calculates momentum and kinetic energy.
    // Requires that one of E, Ek, P be non-zero (only one).
    BDSParticleDefinition* particleDefinition = nullptr;
    BDSIonDefinition* ionDef = nullptr;

    // PDG for ions = 10LZZZAAAI
    if (pdgID > 1000000000) // is an ion
    {   

        G4IonTable* ionTable = particleTable->GetIonTable();
        particleDefGeant = ionTable->GetIon(pdgID);
        if (!particleDefGeant)
        {throw BDSException("BDSXtrackInterface> Ion \"" + std::to_string(pdgID) + "\" not found");}

        G4int ionCharge = ionChargeIn==0 ? particleDefGeant->GetAtomicNumber() : (G4int) ionChargeIn; 
        ionDef = new BDSIonDefinition(particleDefGeant->GetAtomicMass(),
                                      particleDefGeant->GetAtomicNumber(),
                                      ionCharge);

        // correct the particle mass in the case of partially-stripped ions
        G4double mass   = ionTable->GetIonMass(ionDef->Z(), ionDef->A());
        G4double charge = ionDef->Charge(); // correct even if overridden
        mass += ionDef->NElectrons()*G4Electron::Definition()->GetPDGMass();
        // The constructor with a custom name requires the bdsim name of the ion
        std::string bdsimPartName = "ion " + std::to_string(static_cast<int>(ionDef->A())) 
                                  + " " + std::to_string(static_cast<int>(ionDef->Z())) 
                                  + " " + std::to_string(static_cast<int>(charge));

        particleDefinition = new BDSParticleDefinition(bdsimPartName, mass, charge, 0, 
                                                       kineticEnergyIn, momentumIn, 1, ionDef, pdgID);
    }
    else
    {
        particleDefGeant = particleTable->FindParticle(pdgID);
        if (!particleDefGeant)
        {throw BDSException("BDSXtrackInterface> Particle \"" + std::to_string(pdgID) + "\" not found");}
        particleDefinition = new BDSParticleDefinition(particleDefGeant, 0, kineticEnergyIn, momentumIn, 1, nullptr);
    }

    return particleDefinition;
}

BDSParticleDefinition* BDSLinkTrackerInterface::GetReferenceParticleDefinition() const { return referenceParticleDefinition; }



PYBIND11_MODULE(bdslinktrackerinterface, m) {
  py::class_<BDSLinkTrackerInterface>(m,"BDSLinkTrackerInterface")
      .def_static("GetInstance", [](std::string bdsimConfigFile,
                                    int referenceParticlePDG,
                                    double referenceKineticEnergy,
                                    double relativeEnergyCut,
                                    int seed,
                                    int referenceIonCharge,
                                    bool batchMode) {
        auto* obj = BDSLinkTrackerInterface::GetInstance(bdsimConfigFile,
                                                    referenceParticlePDG,
                                                    referenceKineticEnergy,
                                                    relativeEnergyCut,
                                                    seed,
                                                    referenceIonCharge,
                                                    batchMode);
        obj->SetReferenceParticleDefinition(PrepareBDSParticleDefinition(referenceParticlePDG,
                                                                     0,
                                                                     referenceKineticEnergy,
                                                                     referenceIonCharge));
        return obj;
    },
                                    py::arg("bdsimConfigFileIn") = "trackerInterface.gmad",
                                    py::arg("referenceParticlePDG") =11,
                                    py::arg("referenceKineticEnergy") = 100,
                                    py::arg("relativeEnergyCut") = 0.01,
                                    py::arg("seed") = 1234,
                                    py::arg("referenceIonCharge") = 1,
                                    py::arg("batchMode") = true,
                                    py::return_value_policy::reference)
      .def_static("GetInstance", []() {return BDSLinkTrackerInterface::GetInstance();},
                  py::return_value_policy::reference)
      .def("GetBunchLink",&BDSLinkTrackerInterface::GetBunchLink)
      .def("GetBDSIMLink",&BDSLinkTrackerInterface::GetBDSIMLink)
      .def("TrackXSuite",[](BDSLinkTrackerInterface *tracker_interface,
                            int iElement,
                            std::string elementName,
                            py::object particles) {
        py::print("Element BDSIM:",iElement);
        py::print("Particles:", particles);

        py::array_t<double> x = py::cast<py::array_t<double>>(particles.attr("x"));
        py::array_t<double> y = py::cast<py::array_t<double>>(particles.attr("y"));
        py::array_t<double> px = py::cast<py::array_t<double>>(particles.attr("px"));
        py::array_t<double> py = py::cast<py::array_t<double>>(particles.attr("py"));
        py::array_t<double> ct = py::cast<py::array_t<double>>(particles.attr("zeta"));
        py::array_t<double> deltap = py::cast<py::array_t<double>>(particles.attr("delta"));
        py::array_t<double> chi = py::cast<py::array_t<double>>(particles.attr("chi"));
        py::array_t<double> charge_ratio = py::cast<py::array_t<double>>(particles.attr("charge_ratio"));
        py::array_t<double> s = py::cast<py::array_t<double>>(particles.attr("s"));
        py::array_t<int64_t> pdgid = py::cast<py::array_t<int64_t>>(particles.attr("pdg_id"));
        py::array_t<int64_t> trackid = py::cast<py::array_t<int64_t>>(particles.attr("particle_id"));
        py::array_t<int64_t> state = py::cast<py::array_t<int64_t>>(particles.attr("state"));
        py::array_t<int64_t> at_element = py::cast<py::array_t<int64_t>>(particles.attr("at_element"));
        py::array_t<int64_t> at_turn = py::cast<py::array_t<int64_t>>(particles.attr("at_turn"));

        // fill BDSLinkBunch
        auto bunch = tracker_interface->GetBunchLink();
        auto ref = tracker_interface->GetReferenceParticleDefinition();
        // bunch->AddParticle();

        // Fast pointers to numpy arrays
        auto get = [&](const char* name) {
            return py::cast<py::array_t<double>>(particles.attr(name)).mutable_data();
        };
        auto geti = [&](const char* name) {
            return py::cast<py::array_t<int64_t>>(particles.attr(name)).mutable_data();
        };

        double* x_ptr = get("x");
        double* y_ptr = get("y");
        double* px_ptr = get("px");
        double* py_ptr = get("py");
        double* zeta_ptr = get("zeta");
        double* delta_ptr = get("delta");
        double* chi_ptr = get("chi");
        double* charge_ratio_ptr = get("charge_ratio");
        double* s_ptr = get("s");
        int64_t* pdgid_ptr = geti("pdg_id");
        int64_t* trackid_ptr = geti("particle_id");
        int64_t* state_ptr = geti("state");
        int64_t* at_element_ptr = geti("at_element");
        int64_t* at_turn_ptr = geti("at_turn");

//#ifdef CACHE_PTRS
//        double *x_ptr = make_ptr<double>(x);
//        double *y_ptr = make_ptr<double>(y);
//        double *px_ptr = make_ptr<double>(px);
//        double *py_ptr = make_ptr<double>(py);
//        double *zeta_ptr = make_ptr<double>(zeta);
//        double *delta_ptr = make_ptr<double>(delta);
//        double *chi_ptr = make_ptr<double>(chi);
//        double *charge_ratio_ptr = make_ptr<double>(charge_ratio);
//        double *s_ptr = make_ptr<double>(s);
//        double *pdg_id_ptr = make_ptr<double>(pdg_id);
//        double *particle_id_ptr = make_ptr<double>(particle_id);
//        double *state_ptr = make_ptr<double>(state);
//        double *at_element_ptr = make_ptr<double>(at_element);
//        double *at_turn_ptr = make_ptr<double>(at_turn);
//#endif
//
//
//        double* x_ptr = x.mutable_data();
//        double* y_ptr = y.mutable_data();
//        double* px_ptr = px.mutable_data();
//        double* py_ptr = py.mutable_data();
//        double* ct_ptr = ct.mutable_data();
//        double* deltap_ptr = deltap.mutable_data();
//        double* chi_ptr = chi.mutable_data();
//        double* charge_ratio_ptr = charge_ratio.mutable_data();
//        double* s_ptr = s.mutable_data();
//        int64_t* pdgid_ptr = pdgid.mutable_data();
//        int64_t* trackid_ptr = trackid.mutable_data();
//        int64_t* state_ptr = state.mutable_data();
//
//        auto bunch = tracker_interface->GetBunchLink();
//        auto refParticleDefinition = tracker_interface->GetReferenceParticle();  // Or however you access it

        ssize_t n = py::cast<py::array_t<double>>(particles.attr("x")).size();
        auto& active_state = tracker_interface->GetParticleActiveState();

        for (ssize_t i = 0; i < n; ++i) {
            if (state_ptr[i] < 1) {  // there can be different types of alive states in xsuite
                active_state.push_back(false);
                continue;
            }
            else {
                G4double q = charge_ratio_ptr[i] * ref->Charge();
                G4double mass_ratio = charge_ratio_ptr[i] / chi_ptr[i];
                G4double p = ref->Momentum() * (delta_ptr[i] + 1) * mass_ratio;
        
                auto partDef = PrepareBDSParticleDefinition(pdgid_ptr[i], p, 0, q);
                G4double t = - zeta_ptr[i] * CLHEP::m / (ref->Beta() * CLHEP::c_light);
        
                G4double oneplusdelta = (1 + delta_ptr[i]) * mass_ratio;
                G4double xp = px_ptr[i] / oneplusdelta;
                G4double yp = py_ptr[i] / oneplusdelta;
                G4double zp = BDSBunch::CalculateZp(xp, yp, 1);
        
                BDSParticleCoordsFull coords(x_ptr[i] * CLHEP::m,
                                             y_ptr[i] * CLHEP::m,
                                             0,
                                             xp, yp, zp,
                                             t,
                                             0,
                                             partDef->TotalEnergy(),
                                             1);
                active_state.push_back(true);        
                bunch->AddParticle(partDef, coords, trackid_ptr[i], trackid_ptr[i]);
            }
        }


        // run n particles
        auto link = tracker_interface->GetBDSIMLink();
        std::cout << "beam on with: " << bunch->Size() << std::endl;
        selectCollimator(link, elementName);
        link->BeamOn((G4int)bunch->Size());
        std::cout << "finished tracking" << std::endl;
        // get sampler data
        const BDSHitsCollectionSamplerLink* hits = link->SamplerHits();
        std::cout << "sampled hits: " << hits->GetSize() << std::endl;
        size_t hitsCount = 0;
        if (hits)
        {
            hitsCount = hits->GetSize();
        }
        else
        {
            // There were no hits - check if there were any active particles at all coming in
            if (!bunch->Size())
            {
                bunch->GetNextParticleLocal();
            }
        }
        std::cout << "hits: " << hits << std::endl;
        std::cout << "hitsCount: " << hitsCount << std::endl;
        // Count the number of secondary particles
        size_t secondaryCount = 0;
        for (size_t i = 0; i < hitsCount; i++)
        {
            auto hit = (*hits)[i];
            if (hit->externalParticleID != hit->externalParentID) { secondaryCount++; }
        }

        // The output arrays have slots for all particles, regardless of lost or not, and for secondary particles
        size_t output_size = x.size();
        
    auto x_out = py::array(py::buffer_info(
            nullptr,            /* Pointer to data (nullptr -> ask NumPy to allocate!) */
            sizeof(double),     /* Size of one item */
            py::format_descriptor<double>::value, /* Buffer format */
            1,          /* How many dimensions? */
            { output_size },  /* Number of elements for each dimension */
            { sizeof(double) }  /* Strides for each dimension */
    ));

    auto y_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto px_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto py_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto zeta_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto delta_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));         

    auto s_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto pdgid_out = py::array(py::buffer_info(
            nullptr,
            sizeof(int64_t),
            py::format_descriptor<int64_t>::value,
            1,
            { output_size },
            { sizeof(int64_t) }
    ));

    auto trackid_out = py::array(py::buffer_info(
            nullptr,
            sizeof(int64_t),
            py::format_descriptor<int64_t>::value,
            1,
            { output_size },
            { sizeof(int64_t) }
    ));

    auto state_out = py::array(py::buffer_info(
            nullptr,
            sizeof(int64_t),
            py::format_descriptor<int64_t>::value,
            1,
            { output_size },
            { sizeof(int64_t) }
    ));

    auto at_element_out = py::array(py::buffer_info(
            nullptr,
            sizeof(int64_t),
            py::format_descriptor<int64_t>::value,
            1,
            { output_size },
            { sizeof(int64_t) }
    ));

    auto at_turn_out = py::array(py::buffer_info(
            nullptr,
            sizeof(int64_t),
            py::format_descriptor<int64_t>::value,
            1,
            { output_size },
            { sizeof(int64_t) }
    ));      

    auto massratio_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto chargeratio_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));  

    auto x_prod_buf = x_out.request();
    auto y_prod_buf = y_out.request();
    auto px_prod_buf = px_out.request();
    auto py_prod_buf = py_out.request();
    auto zeta_prod_buf = zeta_out.request();
    auto delta_prod_buf = delta_out.request();
    auto s_prod_buf = s_out.request();
    auto pdgid_prod_buf = pdgid_out.request();
    auto trackid_prod_buf = trackid_out.request();
    auto state_prod_buf = state_out.request();
    auto at_element_prod_buf = at_element_out.request();
    auto at_turn_prod_buf = at_turn_out.request();
    auto massratio_prod_buf = massratio_out.request();
    auto chargeratio_prod_buf = chargeratio_out.request();

    auto *x_prod_ptr = (double *) x_prod_buf.ptr;
    auto *y_prod_ptr = (double *) y_prod_buf.ptr;
    auto *px_prod_ptr = (double *) px_prod_buf.ptr;
    auto *py_prod_ptr = (double *) py_prod_buf.ptr;
    auto *zeta_prod_ptr = (double *) zeta_prod_buf.ptr;
    auto *delta_prod_ptr = (double *) delta_prod_buf.ptr;
    auto *s_prod_ptr = (double *) s_prod_buf.ptr;
    auto *pdgid_prod_ptr = (int64_t *) pdgid_prod_buf.ptr;
    auto *trackid_prod_ptr = (int64_t *) trackid_prod_buf.ptr;
    auto *state_prod_ptr = (int64_t *) state_prod_buf.ptr;
    auto *at_element_prod_ptr = (int64_t *) at_element_prod_buf.ptr;
    auto *at_turn_prod_ptr = (int64_t *) at_turn_prod_buf.ptr;
    auto *massratio_prod_ptr = (double *) massratio_prod_buf.ptr;
    auto *chargeratio_prod_ptr = (double *) chargeratio_prod_buf.ptr;

    for (int i=0; i<n; i++) {
        x_prod_ptr[i] = x_ptr[i];
        y_prod_ptr[i] = y_ptr[i];
        px_prod_ptr[i] = px_ptr[i];
        py_prod_ptr[i] = py_ptr[i];
        zeta_prod_ptr[i] = zeta_ptr[i];
        delta_prod_ptr[i] = delta_ptr[i];
        s_prod_ptr[i] = s_ptr[i];
        pdgid_prod_ptr[i] = pdgid_ptr[i];
        trackid_prod_ptr[i] = trackid_ptr[i];
        state_prod_ptr[i] = state_ptr[i];

        //at_element_prod_ptr[i] = at_element_ptr[i];
        //at_turn_prod_ptr[i] = at_turn_ptr[i];
        massratio_prod_ptr[i] = charge_ratio_ptr[i]/chi_ptr[i];
        chargeratio_prod_ptr[i] = charge_ratio_ptr[i];
    }

    // Loop through the particles in the *original* bunch - the primaries
    size_t hits_index = 0;
    bool prim_survied = false;
    // double sum_deltaplusone_sec = 0.0;
    double sum_secondary_energy = 0.0;

    std::cout << "before loop" << std::endl;
    size_t prod_write_index = active_state.size();
    for (size_t i=0; i < active_state.size(); i++){
        std::cout << "start of loop" << std::endl;

        if (!active_state.at(i)){
            continue; // This was an inactive particle that hasn't been processed, do not change it
        }
        std::cout << "loop1" << std::endl;

        //auto xtrack_part = particles[i]; // Get the cached coordinates of the original xtrack particle

        auto part = bunch->GetNextParticle(); // Advance through the bunch
        auto prim_part_id = bunch->CurrentExternalParticleID(); // Get the ID of the primary particle

        std::cout << "part coords: " << part << std::endl; 
        // Now start looping over the hits - the particles to be returned to the tracker
        // These can be primary or secondary particles. Each primary can produce 0, 1, or 2+ products
        // The products need to be sorted to keep the array order - surviving primary particles are all
        // filled in first. If a primary didn't survive, keep the original coordinates and make it inactive.
        // The hits are ordered by primary event, so just need one loop.
        while (hits_index < hitsCount)
        {
            std::cout << "loop2" << std::endl;
            BDSHitSamplerLink* hit = (*hits)[hits_index];
            if (hit->externalParentID != prim_part_id) { // The hits corresponding to the current primary are exhausted
                break;
            }
            std::cout << "loop3" << std::endl;

            const BDSParticleCoordsFull &coords = hit->coords;

            double qratio = hit->charge / ref->Charge();
            double mratio = hit->mass / ref->Mass();

            double dp = (hit->momentum / mratio - ref->Momentum()) / ref->Momentum();

            double collLength = link->GetArcLengthOfLinkElement(iElement);
            /// Need to compensate for the geometry construction in BDSIM
            /// There is a safety margin that is added to the collimator legnth
            double collMargin = 2.5 * BDSSamplerCustom::ChordLength();
            double zt = ref->Beta() * CLHEP::c_light * ((collLength + collMargin) / (CLHEP::c_light * ref->Beta()) - coords.T);

            double oneplusdelta = (1 + dp) * mratio;

            auto track_id = hit->externalParticleID;
            auto parent_id = hit->externalParentID;
            auto pdg_id = hit->pdgID;

            std::cout << "track_id: " << track_id << ", parent_id: " << parent_id << ", pdg_id: " << pdg_id << std::endl;
            if (track_id == parent_id){
                std::cout << "loop4" << std::endl;
                // This is a primary particle as its parent is itself
                prim_survied = true;

                x_prod_ptr[i] = coords.x / CLHEP::m;
                y_prod_ptr[i] = coords.y / CLHEP::m;
                px_prod_ptr[i] = coords.xp * oneplusdelta; // convert back to px proper
                py_prod_ptr[i] = coords.yp * oneplusdelta;
                zeta_prod_ptr[i] = zt / CLHEP::m;
                delta_prod_ptr[i] = dp;
                s_prod_ptr[i] = s_ptr[i] + collLength / CLHEP::m;
                pdgid_prod_ptr[i] = pdg_id;
                //trackid_ptr[i] = track_id; // Don't touch the primary particle id
                trackid_prod_ptr[i] = track_id;
                state_prod_ptr[i] = 1; // active

                at_element_prod_ptr[i] = at_element_ptr[i]; // active
                at_turn_prod_ptr[i] = at_turn_ptr[i]; // active
                massratio_prod_ptr[i] = mratio;
                chargeratio_prod_ptr[i] = qratio;
            }
            else
            {
                std::cout << "loop5" << std::endl;
                // Secondary particles are populated in newly allocated arrays
                x_prod_ptr[prod_write_index] = coords.x / CLHEP::m;
                y_prod_ptr[prod_write_index] = coords.y / CLHEP::m;
                px_prod_ptr[prod_write_index] = coords.xp * oneplusdelta; // convert back to px proper
                py_prod_ptr[prod_write_index] = coords.yp * oneplusdelta; // convert back to py proper;
                zeta_prod_ptr[prod_write_index] = zt / CLHEP::m;
                delta_prod_ptr[prod_write_index] = dp;
                s_prod_ptr[prod_write_index] = s_ptr[i] + collLength / CLHEP::m;
                pdgid_prod_ptr[prod_write_index] = pdg_id;
                trackid_prod_ptr[prod_write_index] = parent_id;
                state_prod_ptr[prod_write_index] = 1; // active
                at_element_prod_ptr[prod_write_index] = at_element_ptr[i]; // active
                at_turn_prod_ptr[prod_write_index] = at_turn_ptr[i]; // active
                massratio_prod_ptr[prod_write_index] = mratio;
                chargeratio_prod_ptr[prod_write_index] = qratio;

                sum_secondary_energy += std::sqrt(std::pow(hit->momentum,2) + std::pow(hit->mass,2));
                prod_write_index++;
            }

            hits_index++;
        }

        if (!prim_survied) // Primary didn't survive - set inactive
        {
            std::cout << "loop6" << std::endl;
            std::cout << "i = " << i << "/" << n << ", state = " << state_prod_ptr[i] << std::endl;
            state_prod_ptr[i] = -333; // inactive
            std::cout << "i = " << i << "/" << n << ", state = " << state_prod_ptr[i] << std::endl;
            // Correct the energy of the lost primary particle to account for the production of secondaries
            // The effective delta is such that the lost particle has the effective delta
            // which corresponds to the energy in - energy out for this primary

            // reconstruct the incoming primary particle energy
            G4double qprim = charge_ratio_ptr[i] * ref->Charge();
            G4double mass_ratio_prim = charge_ratio_ptr[i] / chi_ptr[i];
            G4double p_prim = ref->Momentum() * (delta_ptr[i] + 1) * mass_ratio_prim;
            G4double mass_prim = mass_ratio_prim * ref->Mass();
            G4double energy_prim = std::sqrt(std::pow(p_prim, 2) + std::pow(mass_prim, 2));

            // compute the effective delta
            G4double energy_diff = energy_prim - sum_secondary_energy;
            G4double p_eff;
            G4double squared_energy_mass_diff = std::pow(energy_diff, 2) - std::pow(mass_prim, 2);
            if (squared_energy_mass_diff < 0){
                // This means that the total energy escaping includes part of the rest mass of 
                // the primary. Tolerate the error for now, as otherwise need to adjust also the
                // mass and PDG id of the lost primary particle
                p_eff = 0;
            }
            else
            {
                p_eff = std::sqrt(squared_energy_mass_diff);
            }
            G4double delta_eff = (p_eff / mass_ratio_prim - ref->Momentum()) / ref->Momentum();
            delta_prod_ptr[i] = delta_eff;
        }
        prim_survied = false; // reset for next particle
        sum_secondary_energy = 0.0;
    }

    auto result = py::dict();

    result["s"] = s_out;
    result["x"] = x_out;
    result["px"] = px_out;
    result["y"] = y_out;
    result["py"] = py_out;
    result["zeta"] = zeta_out;
    result["delta"] = delta_out;
    result["pdg_id"] = pdgid_out;
    result["at_element"] = at_element_out;
    result["at_turn"] = at_turn_out;
    result["mass_ratio"] = massratio_out;
    result["charge_ratio"] = chargeratio_out;
    result["parent_particle_id"] = trackid_out;
    result["state"] = state_out;

    std::cout << "result: " << result << std::endl;

        // change state of existing particles
        set_element<double>(x, 0, 1000);      
        // add extra particles of products

    auto state_view = state.mutable_unchecked<1>();
    for (int i=0; i<n; i++) {
        std::cout << "i = " << i << "/" << n << ", x = " << x_prod_ptr[i] << std::endl;
        set_element<double>(s, i, s_prod_ptr[i]);
        set_element<double>(x, i, x_prod_ptr[i]);
        set_element<double>(px, i, px_prod_ptr[i]);
        set_element<double>(y, i, y_prod_ptr[i]);
        set_element<double>(py, i, py_prod_ptr[i]);
        set_element<double>(ct, i, zeta_prod_ptr[i]);
        set_element<double>(deltap, i, delta_prod_ptr[i]);
        set_element<int64_t>(pdgid, i, pdgid_prod_ptr[i]);
        set_element<int64_t>(at_element, i, at_element_prod_ptr[i]);
        set_element<int64_t>(at_turn, i, at_turn_prod_ptr[i]);
//        set_element<double>(mass_ratio, i, massratio_prod_ptr[i]);
        set_element<double>(charge_ratio, i, chargeratio_prod_ptr[i]);
        set_element<int64_t>(trackid, i, trackid_prod_ptr[i]);
        set_element<int64_t>(state, i, state_prod_ptr[i]);
        //state_view(i) = state_prod_ptr[i];
    }        
        set_element<double>(x, 1, 100);      
        set_element<double>(x, 2, 10);
        set_element<double>(s, 0, -1000);      
        set_element<double>(x, 0, -1000);      
        set_element<double>(px, 0, -1000);      
        set_element<double>(y, 0, -1000);      
        set_element<double>(py, 0, -1000);      
        set_element<double>(ct, 0, -1000);      
        set_element<double>(deltap, 0, -1000);      
        set_element<int64_t>(pdgid, 0, -1000);      
        set_element<int64_t>(at_element, 0, -1000);      
        set_element<int64_t>(at_turn, 0, -1000);      
        set_element<double>(charge_ratio, 0, -1000);      
        set_element<int64_t>(trackid, 0, -1000);      
        set_element<int64_t>(state, 0, -1000);      

    return result;

        // clean BDSLinkBunch of particles
      });
}
