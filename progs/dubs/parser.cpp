#include "parser.hpp"

#ifndef _MSVC_LANG
#include <sys/stat.h>
static void mkdir(const std::string& dir) {
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}
#else
#include <direct.h>
static void mkdir(const std::string& dir) {
    _mkdir(dir.c_str());
}
#endif

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

    std::string to_string() const {
        return std::string(begin(), end());
    }

  private:
    const char* begin_;
    size_t length_;
};

static bool is_whitespace(char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

static std::vector<string_view> split_string(string_view string) {
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
static string_view trim(string_view string) {
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

static const std::map<std::string, DUBSParser::TAG_TYPE> TAGS = {
    {"@<no_align_sm_ligands>", DUBSParser::TAG_TYPE::LIGAND_NO_ALIGN},
    {"@<reference>", DUBSParser::TAG_TYPE::REFERENCE},
    {"@<align_prot>", DUBSParser::TAG_TYPE::PROTEIN_ALIGN},
    {"@<align_sm_ligands>", DUBSParser::TAG_TYPE::LIGAND_ALIGN},
    {"@<align_non_sm_ligands>", DUBSParser::TAG_TYPE::PEPTIDE_ALIGN},
    {"@<end>", DUBSParser::TAG_TYPE::END}};

const lemon::ResidueNameSet& DUBSParser::ligands(const std::string& entry) const {
    auto rns = entries_to_rns_.find(entry);
    if (rns == entries_to_rns_.end()) {
        return blank_rns_;
    }

    return rns->second;
}

void DUBSParser::make_directories(const std::string& output_dir) const {
    for (auto& kv : reference_to_name_) {
        if (kv.second.empty()) {
            continue;
        }
        mkdir(output_dir + kv.second);
    }
}

std::string DUBSParser::dump() const {

    std::ostringstream oss;

    oss << "@<no_align_sm_ligands>\n";
    for (auto& current_entry : entries_to_tag_) {
        if (current_entry.second != TAG_TYPE::LIGAND_NO_ALIGN) {
            continue;
        }

        oss << current_entry.first;

        for (auto& ligand : ligands(current_entry.first)) {
            oss << " " << ligand;
        }

        oss << "\n";
    }

    oss << "@<end>\n\n";

    for (auto& ref : reference_entries_) {
        auto ref_name = name(ref.first);
        if (!ref_name.empty()) {
            oss << std::string(ref_name.begin(), ref_name.end()) << "\n";
        }

        oss << "@<reference>\n" << ref.first << " ";
        auto ref_path = path(ref.first);
        oss << std::string(ref_path.begin(), ref_path.end()) << "\n";

        oss << "@<align_prot>\n";
        for (auto& current_entry : entries_to_tag_) {
            if (current_entry.second != TAG_TYPE::PROTEIN_ALIGN) {
                continue;
            }

            if (reference(current_entry.first) != ref.first) {
                continue;
            }

            oss << current_entry.first << "\n";
        }

        oss << "@<align_sm_ligands>\n";
        for (auto& current_entry : entries_to_tag_) {
            if (current_entry.second != TAG_TYPE::LIGAND_ALIGN) {
                continue;
            }

            if (reference(current_entry.first) != ref.first) {
                continue;
            }

            oss << current_entry.first;

            for (auto& ligand : ligands(current_entry.first)) {
                oss << " " << ligand;
            }

            oss << "\n";
        }

        oss << "@<end>\n\n";
    }

    return oss.str();
}


const std::string& DUBSParser::get_mapping(const std::map<std::string, std::string>& map,
                        const std::string& e) const {
    auto needle = map.find(e);
    if (needle == map.end()) {
        return blank_;
    }

    return needle->second;
}

static DUBSParser::TAG_TYPE get_tag(const std::string& s) {
    auto trimmed = trim(s);
        
    auto a = std::string(trimmed.begin(), trimmed.end());

    // strip spaces
    std::transform(a.begin(), a.end(), a.begin(), ::tolower);

    auto tag = TAGS.find(a);
    if (tag == TAGS.end()) {
        return DUBSParser::TAG_TYPE::NONE;
    }

    return tag->second;
}

void DUBSParser::parse_stream(std::istream& i) {
    std::string line;

    while (std::getline(i, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        auto tag = get_tag(line);
        switch (tag) {
        case TAG_TYPE::NONE:
            break;
        case TAG_TYPE::REFERENCE:
            std::getline(i, line);
            parse_reference(std::move(line));
            current_tag_ = TAG_TYPE::NONE;
            continue;
        case TAG_TYPE::END:
            last_line_ = "";
            current_reference_ = "";
            current_tag_ = TAG_TYPE::NONE;
            continue;
        case TAG_TYPE::PEPTIDE_ALIGN:
        case TAG_TYPE::PROTEIN_ALIGN:
        case TAG_TYPE::LIGAND_ALIGN:
        case TAG_TYPE::LIGAND_NO_ALIGN:
        default: // <- Shouldn't be needed
            current_tag_ = tag;
            continue;
            break;
        }

        if (current_tag_ == TAG_TYPE::NONE) {
            last_line_ = std::move(line);
            continue;
        }

        if (current_tag_ == TAG_TYPE::PEPTIDE_ALIGN ||
            current_tag_ == TAG_TYPE::PROTEIN_ALIGN ||
            current_tag_ == TAG_TYPE::LIGAND_ALIGN ||
            current_tag_ == TAG_TYPE::LIGAND_NO_ALIGN) {
                
            parse_complex(std::move(line));
        }
    }
}

void DUBSParser::parse_reference(std::string line) {
    auto split = split_string(trim(line));

    current_reference_ = split[0].to_string();

    entries_to_use_.insert(current_reference_);

    std::transform(current_reference_.begin(), current_reference_.end(),
                    current_reference_.begin(), ::toupper);

    auto reference_path = split[1].to_string();
    reference_to_path_[current_reference_] = reference_path;

    entries_to_tag_[current_reference_] = TAG_TYPE::REFERENCE;

    if (!last_line_.empty()) {
        reference_to_name_[current_reference_] = last_line_;
    }

    chemfiles::Trajectory trj(reference_path, 'r');
    reference_to_structure_[current_reference_] = std::move(trj.read());

    if (split.size() >= 3) {
        auto ligand = split[2].to_string();
        entries_to_rns_[current_reference_] = lemon::ResidueNameSet({ligand});
    }
}

void DUBSParser::parse_complex(std::string line) {
    auto split = split_string(trim(line));
    auto entry = split[0].to_string();

    std::transform(entry.begin(), entry.end(), entry.begin(), ::toupper);

    entries_to_use_.insert(entry);
    entries_to_tag_[entry] = current_tag_;

    if (current_tag_ != TAG_TYPE::LIGAND_NO_ALIGN) {

        if (current_reference_.empty()) {
            // throw an error
        }

        entries_to_reference_[entry] = current_reference_;
        reference_entries_[current_reference_].push_back(entry);
    }

    if (current_tag_ == TAG_TYPE::LIGAND_ALIGN ||
        current_tag_ == TAG_TYPE::LIGAND_NO_ALIGN ) {
        auto ligand = split[1].to_string();
        entries_to_rns_[entry] = lemon::ResidueNameSet({ligand});
        // TODO Suppport multiple ligands
    }
}
