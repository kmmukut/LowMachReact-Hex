#include <iostream>
#include <string>
#include <vector>
#include <memory>

// Cantera headers
#include "cantera/core.h"
#include "cantera/transport.h"
#include <algorithm>

// Global pointer to the Cantera solution
static std::shared_ptr<Cantera::Solution> sol;
static std::shared_ptr<Cantera::ThermoPhase> gas;
static std::shared_ptr<Cantera::Transport> trans;

extern "C" {

    // Get number of species in mechanism
    int cantera_get_species_count_c() {
        if (!gas) return 0;
        return (int)gas->nSpecies();
    }

    // Get name of species by index
    void cantera_get_species_name_c(int k, char* name_out, int name_len) {
        if (!gas || k < 0 || k >= (int)gas->nSpecies()) return;
        std::string name = gas->speciesName(k);
        // Copy and pad with spaces for Fortran
        int len = std::min((int)name.length(), name_len);
        for (int i = 0; i < len; ++i) name_out[i] = name[i];
        for (int i = len; i < name_len; ++i) name_out[i] = ' ';
    }

    // Initialize the Cantera mechanism
    void cantera_init_c(const char* mech_file, int nspecies, const char* species_names_flat, int name_len) {
        try {
            std::string mech_str(mech_file);
            // Trim trailing spaces from mech_file
            mech_str.erase(mech_str.find_last_not_of(" ") + 1);

            std::cout << "Cantera Bridge: Initializing with mechanism: " << mech_str << std::endl;

            // Load the gas mixture
            sol = Cantera::newSolution(mech_str);
            gas = sol->thermo();
            
            // Create transport manager
            trans = sol->transport();

            std::cout << "Cantera Bridge: Initialization successful." << std::endl;
            
        } catch (Cantera::CanteraError& err) {
            std::cerr << "Cantera Error: " << err.what() << std::endl;
            exit(1);
        } catch (std::exception& e) {
            std::cerr << "Standard Exception: " << e.what() << std::endl;
            exit(1);
        }
    }

    // Compute transport properties
    void cantera_update_transport_c(int ncells, double* T, double* P, int nspecies, double* Y_in, 
                                    double* mu_out, double* diff_out, 
                                    const char* species_names_flat, int name_len) {
        if (!gas || !trans) {
            std::cerr << "Cantera Bridge: Gas not initialized before update_transport!" << std::endl;
            exit(1);
        }

        try {
            int cantera_nsp = gas->nSpecies();
            std::vector<double> Y_cantera(cantera_nsp, 0.0);
            std::vector<double> diff_cantera(cantera_nsp, 0.0);

            // Get mapping between solver species and Cantera species
            std::vector<int> sp_map(nspecies, -1);
            int bath_gas_index = gas->speciesIndex("N2"); // Default bath gas

            for (int k = 0; k < nspecies; ++k) {
                std::string sp_name(species_names_flat + k * name_len, name_len);
                sp_name.erase(sp_name.find_last_not_of(" ") + 1);
                
                int c_idx = gas->speciesIndex(sp_name);
                if (c_idx != Cantera::npos) {
                    sp_map[k] = c_idx;
                } else {
                    std::cerr << "Cantera Bridge: Species " << sp_name << " not found in mechanism!" << std::endl;
                }
            }

            for (int c = 0; c < ncells; ++c) {
                // Reset Cantera mass fractions
                std::fill(Y_cantera.begin(), Y_cantera.end(), 0.0);
                double sum_Y = 0.0;

                // Set solver species
                for (int k = 0; k < nspecies; ++k) {
                    if (sp_map[k] >= 0) {
                        double y_val = Y_in[c * nspecies + k];
                        // Ensure non-negative mass fractions
                        if (y_val < 0.0) y_val = 0.0;
                        Y_cantera[sp_map[k]] = y_val;
                        sum_Y += y_val;
                    }
                }

                // Put the remaining mass into bath gas
                if (bath_gas_index >= 0) {
                    Y_cantera[bath_gas_index] = std::max(0.0, 1.0 - sum_Y);
                } else {
                    // If no N2, just normalize what we have
                }

                // Set Cantera state
                gas->setState_TPY(T[c], P[c], Y_cantera.data());

                // Get dynamic viscosity
                mu_out[c] = trans->viscosity();

                // Get mixture-averaged diffusion coefficients
                trans->getMixDiffCoeffs(diff_cantera.data());

                // Map back to solver species array
                for (int k = 0; k < nspecies; ++k) {
                    if (sp_map[k] >= 0) {
                        diff_out[c * nspecies + k] = diff_cantera[sp_map[k]];
                    } else {
                        diff_out[c * nspecies + k] = 0.0;
                    }
                }
            }

        } catch (Cantera::CanteraError& err) {
            std::cerr << "Cantera Error: " << err.what() << std::endl;
            exit(1);
        }
    }

} // extern "C"
