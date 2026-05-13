#include <iostream>
#include <string>
#include <vector>
#include <memory>

// Cantera headers
#include <cantera/core.h>
#include <cantera/transport/TransportFactory.h>
#include <algorithm>
#include <cmath>

// Global pointer to the Cantera solution
static std::shared_ptr<Cantera::Solution> sol;
static std::shared_ptr<Cantera::ThermoPhase> gas;
static std::shared_ptr<Cantera::Transport> trans;

// Build a full Cantera mass-fraction vector from the solver species list.
// If the solver list is empty, or if the sum is below one, the remaining
// mass fraction is assigned to N2 when present. This is the fallback for
// non-reacting no-species thermo paths.
static void build_cantera_mass_fractions(
    Cantera::ThermoPhase& gas_phase,
    int nspecies,
    const double* Y_in,
    int cell_index,
    const std::vector<int>& sp_map,
    std::vector<double>& Y_cantera)
{
    std::fill(Y_cantera.begin(), Y_cantera.end(), 0.0);

    int bath_gas_index = (int)gas_phase.speciesIndex("N2");
    double sum_Y = 0.0;

    for (int k = 0; k < nspecies; ++k) {
        if (sp_map[k] >= 0) {
            double y_val = Y_in[cell_index * nspecies + k];
            if (y_val < 0.0) {
                y_val = 0.0;
            }
            Y_cantera[sp_map[k]] = y_val;
            sum_Y += y_val;
        }
    }

    if (bath_gas_index >= 0) {
        Y_cantera[bath_gas_index] += std::max(0.0, 1.0 - sum_Y);
    }

    double total = 0.0;
    for (double y : Y_cantera) {
        total += y;
    }

    if (total <= 0.0) {
        if (bath_gas_index >= 0) {
            Y_cantera[bath_gas_index] = 1.0;
        } else if (!Y_cantera.empty()) {
            Y_cantera[0] = 1.0;
        }
    }
}

// Map the solver species names to Cantera species indices.
static std::vector<int> build_species_map(
    Cantera::ThermoPhase& gas_phase,
    int nspecies,
    const char* species_names_flat,
    int name_len)
{
    std::vector<int> sp_map(nspecies, -1);

    for (int k = 0; k < nspecies; ++k) {
        std::string sp_name(species_names_flat + k * name_len, name_len);
        size_t last = sp_name.find_last_not_of(" ");
        if (last == std::string::npos) {
            sp_name = "";
        } else {
            sp_name.erase(last + 1);
        }

        size_t c_idx = gas_phase.speciesIndex(sp_name);
        if (c_idx != Cantera::npos) {
            sp_map[k] = (int)c_idx;
        } else {
            std::cerr << "Cantera Bridge: Species " << sp_name
                      << " not found in mechanism!" << std::endl;
        }
    }

    return sp_map;
}


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
        size_t len = std::min(name.length(), (size_t)name_len);
        for (size_t i = 0; i < len; ++i) name_out[i] = name[i];
        for (size_t i = len; i < (size_t)name_len; ++i) name_out[i] = ' ';
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
            static int cached_ncells = 0;
            static int cached_nspecies = 0;
            static std::vector<double> cached_Y;
            static std::vector<double> cached_T;
            static std::vector<double> cached_P;
            static std::vector<double> cached_mu;
            static std::vector<double> cached_diff;

            if (ncells != cached_ncells || nspecies != cached_nspecies) {
                cached_Y.assign(ncells * std::max(1, nspecies), -1.0);
                cached_T.assign(ncells, -1.0);
                cached_P.assign(ncells, -1.0);
                cached_mu.assign(ncells, 0.0);
                cached_diff.assign(ncells * std::max(1, nspecies), 0.0);
                cached_ncells = ncells;
                cached_nspecies = nspecies;
            }

            int cantera_nsp = (int)gas->nSpecies();
            std::vector<double> Y_cantera(cantera_nsp, 0.0);
            std::vector<double> diff_cantera(cantera_nsp, 0.0);

            // Get mapping between solver species and Cantera species
            std::vector<int> sp_map(nspecies, -1);
            int bath_gas_index = (int)gas->speciesIndex("N2"); // Default bath gas

            for (int k = 0; k < nspecies; ++k) {
                std::string sp_name(species_names_flat + k * name_len, name_len);
                sp_name.erase(sp_name.find_last_not_of(" ") + 1);
                
                size_t c_idx = gas->speciesIndex(sp_name);
                if (c_idx != Cantera::npos) {
                    sp_map[k] = (int)c_idx;
                } else {
                    std::cerr << "Cantera Bridge: Species " << sp_name << " not found in mechanism!" << std::endl;
                }
            }

            for (int c = 0; c < ncells; ++c) {
                bool use_cache = true;
                double max_diff = 0.0;
                for (int k = 0; k < nspecies; ++k) {
                    double diff = std::abs(Y_in[c * nspecies + k] - cached_Y[c * nspecies + k]);
                    if (diff > max_diff) max_diff = diff;
                }
                if (std::abs(T[c] - cached_T[c]) > 1.0e-6 ||
                    std::abs(P[c] - cached_P[c]) > 1.0e-3 ||
                    max_diff > 1.0e-5 || cached_mu[c] <= 0.0) {
                    use_cache = false;
                }

                if (use_cache) {
                    mu_out[c] = cached_mu[c];
                    for (int k = 0; k < nspecies; ++k) {
                        diff_out[c * nspecies + k] = cached_diff[c * nspecies + k];
                    }
                    continue; // Skip Cantera evaluation
                }

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

                cached_mu[c] = mu_out[c];

                // Map back to solver species array
                for (int k = 0; k < nspecies; ++k) {
                    if (sp_map[k] >= 0) {
                        diff_out[c * nspecies + k] = diff_cantera[sp_map[k]];
                    } else {
                        diff_out[c * nspecies + k] = 0.0;
                    }
                    cached_diff[c * nspecies + k] = diff_out[c * nspecies + k];
                    cached_Y[c * nspecies + k] = Y_in[c * nspecies + k];
                }
                cached_T[c] = T[c];
                cached_P[c] = P[c];
            }

        } catch (Cantera::CanteraError& err) {
            std::cerr << "Cantera Error: " << err.what() << std::endl;
            exit(1);
        }
    }


    // Compute thermodynamic properties from T, P, and composition.
    void cantera_update_thermo_c(int ncells, double* T, double* P, int nspecies, double* Y_in,
                                 double* h_out, double* cp_out, double* lambda_out,
                                 double* rho_thermo_out, double T_ref,
                                 const char* species_names_flat, int name_len) {
        if (!gas || !trans) {
            std::cerr << "Cantera Bridge: Gas not initialized before update_thermo!" << std::endl;
            exit(1);
        }

        try {
            int cantera_nsp = (int)gas->nSpecies();
            std::vector<double> Y_cantera(cantera_nsp, 0.0);
            std::vector<int> sp_map = build_species_map(*gas, nspecies, species_names_flat, name_len);

            for (int c = 0; c < ncells; ++c) {
                build_cantera_mass_fractions(*gas, nspecies, Y_in, c, sp_map, Y_cantera);
                gas->setState_TPY(T[c], P[c], Y_cantera.data());

                const double h_abs = gas->enthalpy_mass();
                cp_out[c] = gas->cp_mass();
                rho_thermo_out[c] = gas->density();
                lambda_out[c] = trans->thermalConductivity();

                // Non-reacting energy uses sensible enthalpy, not
                // total formation-including enthalpy. This avoids artificial
                // heat release when composition changes without reactions.
                gas->setState_TPY(T_ref, P[c], Y_cantera.data());
                const double h_ref = gas->enthalpy_mass();
                h_out[c] = h_abs - h_ref;
            }

        } catch (Cantera::CanteraError& err) {
            std::cerr << "Cantera Error: " << err.what() << std::endl;
            exit(1);
        } catch (std::exception& e) {
            std::cerr << "Standard Exception: " << e.what() << std::endl;
            exit(1);
        }
    }


    // Recover temperature from mixture enthalpy, pressure, and composition.
    void cantera_recover_temperature_from_h_c(int ncells, double* h_in, double* P, int nspecies, double* Y_in,
                                              double* T_out, double T_ref,
                                              const char* species_names_flat, int name_len) {
        if (!gas) {
            std::cerr << "Cantera Bridge: Gas not initialized before recover_temperature!" << std::endl;
            exit(1);
        }

        try {
            int cantera_nsp = (int)gas->nSpecies();
            std::vector<double> Y_cantera(cantera_nsp, 0.0);
            std::vector<int> sp_map = build_species_map(*gas, nspecies, species_names_flat, name_len);

            for (int c = 0; c < ncells; ++c) {
                build_cantera_mass_fractions(*gas, nspecies, Y_in, c, sp_map, Y_cantera);

                // h_in is sensible enthalpy relative to T_ref for the same
                // composition. Convert back to Cantera absolute enthalpy
                // before HP inversion.
                gas->setState_TPY(T_ref, P[c], Y_cantera.data());
                const double h_ref = gas->enthalpy_mass();
                const double h_abs_target = h_in[c] + h_ref;

                gas->setMassFractions(Y_cantera.data());
                gas->setState_HP(h_abs_target, P[c]);
                T_out[c] = gas->temperature();
            }

        } catch (Cantera::CanteraError& err) {
            std::cerr << "Cantera Error: " << err.what() << std::endl;
            exit(1);
        } catch (std::exception& e) {
            std::cerr << "Standard Exception: " << e.what() << std::endl;
            exit(1);
        }
    }


    // Recover temperature from sensible enthalpy and refresh thermo properties
    // in one cell loop.  This preserves the transported h_in array and avoids
    // the previous recover_temperature() + update_thermo() double Cantera pass.
    //
    // A conservative cache is included for cells whose h, p, and solver species
    // mass fractions are unchanged to roundoff-level tolerances.
    void cantera_recover_temperature_and_update_thermo_c(
        int ncells, double* h_in, double* P, int nspecies, double* Y_in,
        double* T_out, double* cp_out, double* lambda_out, double* rho_thermo_out,
        double T_ref, const char* species_names_flat, int name_len)
    {
        if (!gas || !trans) {
            std::cerr << "Cantera Bridge: Gas not initialized before recover_temperature_and_update_thermo!"
                      << std::endl;
            exit(1);
        }

        try {
            static int cached_ncells = 0;
            static int cached_nspecies = 0;
            static std::vector<double> cached_h;
            static std::vector<double> cached_P;
            static std::vector<double> cached_Y;
            static std::vector<double> cached_T;
            static std::vector<double> cached_cp;
            static std::vector<double> cached_lambda;
            static std::vector<double> cached_rho;

            if (ncells != cached_ncells || nspecies != cached_nspecies) {
                cached_h.assign(ncells, 1.0e300);
                cached_P.assign(ncells, 1.0e300);
                cached_Y.assign(ncells * std::max(1, nspecies), 1.0e300);
                cached_T.assign(ncells, 0.0);
                cached_cp.assign(ncells, 0.0);
                cached_lambda.assign(ncells, 0.0);
                cached_rho.assign(ncells, 0.0);
                cached_ncells = ncells;
                cached_nspecies = nspecies;
            }

            const double h_rel_tol = 1.0e-12;
            const double p_abs_tol = 1.0e-3;
            const double y_abs_tol = 1.0e-12;

            int cantera_nsp = (int)gas->nSpecies();
            std::vector<double> Y_cantera(cantera_nsp, 0.0);
            std::vector<int> sp_map = build_species_map(*gas, nspecies, species_names_flat, name_len);

            for (int c = 0; c < ncells; ++c) {
                bool use_cache = (cached_cp[c] > 0.0 && cached_lambda[c] >= 0.0 && cached_rho[c] > 0.0);

                const double h_tol = h_rel_tol * std::max(1.0, std::abs(h_in[c]));
                if (std::abs(h_in[c] - cached_h[c]) > h_tol ||
                    std::abs(P[c] - cached_P[c]) > p_abs_tol) {
                    use_cache = false;
                }

                if (use_cache) {
                    for (int k = 0; k < nspecies; ++k) {
                        if (std::abs(Y_in[c * nspecies + k] - cached_Y[c * nspecies + k]) > y_abs_tol) {
                            use_cache = false;
                            break;
                        }
                    }
                }

                if (use_cache) {
                    T_out[c] = cached_T[c];
                    cp_out[c] = cached_cp[c];
                    lambda_out[c] = cached_lambda[c];
                    rho_thermo_out[c] = cached_rho[c];
                    continue;
                }

                build_cantera_mass_fractions(*gas, nspecies, Y_in, c, sp_map, Y_cantera);

                // h_in is sensible enthalpy relative to T_ref for the same
                // composition. Convert back to Cantera absolute enthalpy
                // before HP inversion.
                gas->setState_TPY(T_ref, P[c], Y_cantera.data());
                const double h_ref = gas->enthalpy_mass();
                const double h_abs_target = h_in[c] + h_ref;

                gas->setMassFractions(Y_cantera.data());
                gas->setState_HP(h_abs_target, P[c]);

                T_out[c] = gas->temperature();
                cp_out[c] = gas->cp_mass();
                rho_thermo_out[c] = gas->density();
                lambda_out[c] = trans->thermalConductivity();

                cached_h[c] = h_in[c];
                cached_P[c] = P[c];
                cached_T[c] = T_out[c];
                cached_cp[c] = cp_out[c];
                cached_lambda[c] = lambda_out[c];
                cached_rho[c] = rho_thermo_out[c];

                for (int k = 0; k < nspecies; ++k) {
                    cached_Y[c * nspecies + k] = Y_in[c * nspecies + k];
                }
            }

        } catch (Cantera::CanteraError& err) {
            std::cerr << "Cantera Error: " << err.what() << std::endl;
            exit(1);
        } catch (std::exception& e) {
            std::cerr << "Standard Exception: " << e.what() << std::endl;
            exit(1);
        }
    }

} // extern "C"
