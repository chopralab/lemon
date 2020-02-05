#ifndef DUBS_PARSER_HPP
#define DUBS_PARSER_HPP

#include "lemon/external/gaurd.hpp"

LEMON_EXTERNAL_FILE_PUSH
#include <chemfiles.hpp>
LEMON_EXTERNAL_FILE_POP

#include "lemon/entries.hpp"

#include <string>
#include <map>

class DUBSParser {
  public:

    void parse(std::istream& i) { parse_stream(i); }

    enum class TAG_TYPE {
        NONE = 0,
        LIGAND_NO_ALIGN,
        REFERENCE,
        PROTEIN_ALIGN,
        LIGAND_ALIGN,
        PEPTIDE_ALIGN,
        END
    };

    static const std::map<std::string, TAG_TYPE> TAGS;

    const lemon::Entries& entries() const { return entries_to_use_; }

    //! The reference for a given entry. Returns a blank string if no reference
    //! is found
    const std::string& reference(const std::string& entry) const {
        return get_mapping(entries_to_reference_, entry);
    }

    //! The name for a reference protein. Returns a blank string if no name was
    //! given for the reference
    const std::string& name(const std::string& reference) const {
        return get_mapping(reference_to_name_, reference);
    }

    //! The filesystem path for a reference protein
    const std::string& path(const std::string& reference) const {
        return get_mapping(reference_to_path_, reference);
    }

    //! The TAG used for a given entry
    DUBSParser::TAG_TYPE tag_type(const std::string& entry) const {
        return entries_to_tag_.at(entry);
    }

    //! Retreive the object representing the reference structure of a given
    //! reference protein
    const chemfiles::Frame& reference_structure(const std::string& reference) const {
        return reference_to_structure_.at(reference);
    }

    //! The ligands specified for a reference protein
    const lemon::ResidueNameSet& ligands(const std::string& entry) const;

    void make_directories(const std::string& output_dir) const;

    std::string dump() const;

  private:

    using string_to_string = std::map<std::string, std::string>;
    using string_to_string_vector = std::map<std::string, std::vector<std::string>>;

    TAG_TYPE current_tag_ = TAG_TYPE::NONE;
    const std::string blank_ = "";
    const lemon::ResidueNameSet blank_rns_ = lemon::ResidueNameSet();

    const std::string& get_mapping(const std::map<std::string, std::string>& map,
                            const std::string& e) const;

    //! The last comment line (could be the name of the current target protein)
    std::string last_line_;

    //************************************************************************
    //* Objects related to all entries
    //************************************************************************

    //! All entries that should be loaded by Lemon
    lemon::Entries entries_to_use_;

    //! The residue name sets for a given entry
    std::map<std::string, lemon::ResidueNameSet> entries_to_rns_;

    //! The tag used to add a given entry
    std::map<std::string, TAG_TYPE> entries_to_tag_;

    //************************************************************************
    //* Objects related to the reference protein
    //************************************************************************

    //! The current reference protein
    std::string current_reference_;

    //! The path to the reference files
    string_to_string reference_to_path_;

    //! The structure of the current reference files
    std::map<std::string, chemfiles::Frame> reference_to_structure_;

    //! The name of the target protein
    string_to_string reference_to_name_;

    //! All entries for a given reference protein
    string_to_string_vector reference_entries_;

    //! The reference name for a given entry
    string_to_string entries_to_reference_;

    void parse_stream(std::istream& i);

    void parse_reference(std::string line);

    void parse_complex(std::string line);
};

#endif // DUBS_PARSER_HPP
