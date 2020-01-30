#include "lemon/launch.hpp"
#include "lemon/lemon.hpp"

#include <map>

#include <algorithm>
#include <iostream>
#include <sstream>

struct string_view {

    string_view(const char* c, size_t length) : begin_(c), length_(length) {}
    string_view(const std::string& s) : begin_(&s[0]), length_(s.length()) {}
    string_view(const string_view& o) = default;

    string_view substr(size_t pos, size_t n) {
        if (pos > size()) {
            throw std::out_of_range("string_view::substr()");
        }

        return string_view(begin() + pos, std::min(n, size() - pos));
    }

    char operator[](size_t pos) { return begin_[pos]; }

    const char* begin() const { return begin_; }
    const char* end() const { return begin_ + length_; }

    size_t length() const { return length_; }
    size_t size() const { return length_; }
    bool empty() const { return length_ == 0; }

  private:
    const char* begin_;
    size_t length_;
};

inline bool is_whitespace(char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

inline std::vector<string_view> split_string(string_view string) {
    std::vector<string_view> elems;
    size_t last = 0;
    for (size_t i = 0; i < string.length(); i++) {
        if (is_whitespace(string[i])) {
            if (last != i) {
                // Don't add empty elements
                elems.push_back(string.substr(last, i - last));
            }
            last = i + 1;
        }
    }

    if (last < string.length()) {
        elems.push_back(string.substr(last, string.length() - last));
    }

    return elems;
}

/// Remove whitespaces at the begining and end of `string`
inline string_view trim(string_view string) {
    auto begin = string.begin();
    auto end = string.end();
    while (begin != end && is_whitespace(*begin)) {
        begin++;
    }

    if (begin != end) {
        end--;
        while (end != begin && is_whitespace(*end)) {
            end--;
        }
        end++;
    }

    return string_view(begin, static_cast<size_t>(end - begin));
}

class DUBSParser {
  public:
    DUBSParser(std::istream& i) { parse_stream(i); }

    enum class TAG_TYPE {
        NONE = 0,
        REFERENCE,
        PROTEIN_ALIGN,
        LIGAND_ALIGN,
        PEPTIDE_ALIGN,
        END
    };

    static const std::map<std::string, TAG_TYPE> TAGS;

    const lemon::Entries& entries() const { return entries_to_use_; }

    string_view reference(const std::string& entry) const {
        return get_mapping(entries_to_reference_, entry);
    }

    string_view name(const std::string& reference) const {
        return get_mapping(reference_to_name_, reference);
    }

    string_view path(const std::string& reference) const {
        return get_mapping(reference_to_path_, reference);
    }

    std::string dump() const {

        std::ostringstream oss;

        for (auto& ref : reference_entries_) {
            auto ref_name = name(ref.first);
            if (!ref_name.empty()) {
                oss << std::string(ref_name.begin(), ref_name.end()) << "\n";
            }

            oss << "@<reference>\n" << ref.first << " ";
            auto ref_path = path(ref.first);
            oss << std::string(ref_path.begin(), ref_path.end()) << "\n";

            auto align_prot_iter = reference_align_prot_.find(ref.first);
            if (align_prot_iter != reference_align_prot_.end()) {
                oss << "@<align_prot>\n";
                for (auto& align_prot : align_prot_iter->second) {
                    oss << align_prot << "\n";
                }
            }

            for (auto& entry : ref.second) {
                auto small_molecule_iter = entries_to_rns_.find(entry);
                if (small_molecule_iter != entries_to_rns_.end()) {
                    oss << "@<align_sm_ligands>\n" << ref.first << " ";
                    for (auto& align_sm : small_molecule_iter->second) {
                        oss << align_sm << "\n";
                    }
                }
            }

            oss << "@<end>\n";
        }

        return oss.str();
    }

  private:
    TAG_TYPE current_tag_ = TAG_TYPE::NONE;
    const std::string blank_ = "";

    string_view get_mapping(const std::map<std::string, std::string>& map,
                            const std::string& e) const {
        auto needle = map.find(e);
        if (needle == map.end()) {
            return blank_;
        }

        return needle->second;
    }

    static TAG_TYPE get_tag(const std::string& s) {
        auto trimmed = trim(s);
        
        auto a = std::string(trimmed.begin(), trimmed.end());

        // strip spaces
        std::transform(a.begin(), a.end(), a.begin(), ::tolower);

        auto tag = TAGS.find(a);
        if (tag == TAGS.end()) {
            return TAG_TYPE::NONE;
        }

        return tag->second;
    }

    std::string last_line_;
    std::string current_reference_;
    std::map<std::string, std::string> reference_to_path_;
    std::map<std::string, std::string> reference_to_name_;
    std::map<std::string, std::vector<std::string>> reference_align_prot_;
    std::map<std::string, std::vector<std::string>> reference_entries_;

    lemon::Entries entries_to_use_;
    std::map<std::string, lemon::ResidueNameSet> entries_to_rns_;
    std::map<std::string, std::string> entries_to_reference_;

    void parse_stream(std::istream& i) {
        std::string line;
        while (std::getline(i, line)) {

            if (line == "") {
                continue;
            }

            auto tag = get_tag(line);
            switch (tag) {
            case TAG_TYPE::NONE:
                break;
            case TAG_TYPE::REFERENCE:
                std::getline(i, line);
                parse_reference(std::move(line));
                continue;
            case TAG_TYPE::END:
                last_line_ = "";
                current_reference_ = "";
                current_tag_ = TAG_TYPE::NONE;
                continue;
            case TAG_TYPE::PEPTIDE_ALIGN:
            case TAG_TYPE::PROTEIN_ALIGN:
            case TAG_TYPE::LIGAND_ALIGN:
            default: // <- Shouldn't be needed
                current_tag_ = tag;
                break;
            }

            if (current_tag_ == TAG_TYPE::NONE) {
                last_line_ = std::move(line);
                continue;
            }

            if (current_tag_ == TAG_TYPE::PEPTIDE_ALIGN ||
                current_tag_ == TAG_TYPE::LIGAND_ALIGN ||
                current_tag_ == TAG_TYPE::PROTEIN_ALIGN) {
                auto split = split_string(trim(line));
                auto entry = std::string(split[0].begin(), split[0].end());

                std::transform(entry.begin(), entry.end(), entry.begin(),
                               ::toupper);

                entries_to_use_.insert(entry);

                if (!current_reference_.empty()) {
                    entries_to_reference_[entry] = current_reference_;
                    reference_entries_[current_reference_].push_back(entry);
                }

                if (current_tag_ == TAG_TYPE::PROTEIN_ALIGN) {
                    reference_align_prot_[current_reference_].push_back(entry);
                }

                if (current_tag_ == TAG_TYPE::LIGAND_ALIGN) {
                    auto ligand = std::string(split[1].begin(), split[1].end());
                    entries_to_rns_[entry] = lemon::ResidueNameSet({ligand});
                    // TODO Suppport multiple liands
                }
            }
        }
    }

    void parse_reference(std::string line) {
        auto split = split_string(trim(line));

        current_reference_ = std::string(split[0].begin(), split[0].end());

        std::transform(current_reference_.begin(), current_reference_.end(),
                       current_reference_.begin(), ::toupper);

        auto reference_path = std::string(split[1].begin(), split[1].end());
        reference_to_path_[current_reference_] = reference_path;

        if (!last_line_.empty()) {
            reference_to_name_[current_reference_] = last_line_;
        }
    }
};

const std::map<std::string, DUBSParser::TAG_TYPE> DUBSParser::TAGS = {
    {"<@reference>", TAG_TYPE::REFERENCE},
    {"@<align_prot>", TAG_TYPE::PROTEIN_ALIGN},
    {"@<align_sm_ligands>", TAG_TYPE::LIGAND_ALIGN},
    {"@<align_non_sm_ligands>", TAG_TYPE::PEPTIDE_ALIGN},
    {"@<end>", TAG_TYPE::END}};

int main(int argc, char* argv[]) {
    lemon::Options o;
    std::string outdir = ".";
    auto distance = 6.0;
    o.add_option("--distance,-d", distance,
                 "Largest distance between protein and a small molecule.");
    o.add_option("--outdir,-o", outdir, "output directory");
    o.parse_command_line(argc, argv);

    outdir += "/";

    std::ifstream is(o.entries());

    DUBSParser parser(is);
    std::cout << parser.dump();

    /*
    auto worker = [distance, &outdir](chemfiles::Frame entry,
                                      const std::string& pdbid) {
        // Selection phase
        std::list<size_t> smallm;
        if (lemon::select::specific_residues(entry, smallm, rnms[pdbid]) == 0) {
            return "Skipping " + pdbid;
        }

        // Pruning phase
        lemon::prune::identical_residues(entry, smallm);

        // Output phase
        for (auto resid : smallm) {
            chemfiles::Frame prot;
            chemfiles::Frame lig;
            lemon::separate::protein_and_ligand(entry, resid, distance, prot,
                                                lig);

            auto protfile =
                outdir + pdbid + "_" + lig.get("name")->as_string() + ".pdb";
            auto ligfile =
                outdir + pdbid + "_" + lig.get("name")->as_string() + ".sdf";

            chemfiles::Trajectory prot_traj(protfile, 'w');
            chemfiles::Trajectory lig_traj(ligfile, 'w');
            prot_traj.write(prot);
            lig_traj.write(lig);
            prot_traj.close();
            lig_traj.close();
        }

        return pdbid + std::to_string(smallm.size());
    };

    auto collector = lemon::print_combine(std::cout);
    auto p = o.work_dir();
    auto threads = o.ncpu();

    try {
        lemon::run_parallel(worker, p, collector, threads, entries, entries);
    } catch (std::runtime_error& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }*/

    return 0;
}
