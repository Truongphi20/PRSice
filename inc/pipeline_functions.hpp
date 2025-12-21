#ifndef PIPELINE_FUNCTIONS_HPP
#define PIPELINE_FUNCTIONS_HPP

#include "IITree.h"
#include "genotype.hpp"
#include "genotypefactory.hpp"
#include "misc.hpp"
#include "region.hpp"
#include "reporter.hpp"

#include <exception>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

// ---------- inline functions (OK in header) ----------

inline void print_empty_region(
    const std::string& out,
    const std::vector<std::vector<size_t>>& region_membership,
    const std::vector<std::string>& region_names)
{
    bool has_empty_region = false;
    std::ostringstream empty_regions;
    std::string empty_region_name = out + ".xregion";

    for (size_t region_idx = 2; region_idx < region_membership.size();
         ++region_idx)
    {
        if (region_membership[region_idx].empty())
        {
            has_empty_region = true;
            empty_regions << region_names[region_idx] << '\n';
        }
    }

    if (has_empty_region)
    {
        std::ofstream empty_region_file(empty_region_name);
        if (!empty_region_file)
        {
            throw std::runtime_error(
                "Error: Cannot open file: " + empty_region_name);
        }
        empty_region_file << empty_regions.str();
    }
}

inline void initialize_genotype(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const Commander& commander, Genotype* current_file,
    Reporter& reporter, Genotype* target_file);

inline void initialize_target(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const Commander& commander, Genotype* target_file,
    Reporter& reporter);

inline void initialize_reference(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    const Commander& commander, Genotype* target_file,
    Genotype* reference_file, Reporter& reporter);

inline std::tuple<std::vector<std::string>, size_t>
add_gene_set_info(const Commander& commander,
                  Genotype* target_file,
                  Reporter& reporter);

// ---------- NON-inline declarations ONLY ----------

std::string
print_project_summary(std::vector<size_t>& significant_store);

void print_prsice_header(bool has_prevalence, bool no_regress,
                         std::unique_ptr<std::ostream>& prsice_out);

void print_summary_header(bool has_prevalence, bool run_set_perm,
                          bool run_perm,
                          std::unique_ptr<std::ostream>& summary_file);

#endif // PIPELINE_FUNCTIONS_HPP
