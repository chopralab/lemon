#ifndef LEMON_OPTIONS_HPP
#define LEMON_OPTIONS_HPP

#include "lemon/external/gaurd.hpp"

LEMON_EXTERNAL_FILE_PUSH
// NOLINTNEXTLINE External library uses a define
#define CLI11_EXPERIMENTAL_OPTIONAL 0
#include "lemon/external/CLI11.hpp"
LEMON_EXTERNAL_FILE_POP

namespace lemon {

//! The `Options` class is used to read command line arguments.
//!
//! A majority of **Lemon** workflows begin with the parsing of various command-
//! line arguments. This class helps the user parse the most common arguments
//! such as the directory containing the Hadoop sequence files, the number of
//! cores to use for a given search or a preselected set of entries.
//! Users can add their own options with the `add_options` member function.
class Options : public CLI::App {
  public:
    //! Default constructor to initialize a custom `Options` class
    //!
    //! This constructor is intended for users who wish to add their own custom
    //! options. The *work_dir*, *ncpu*, and *entries* options are added
    //! automatically and additional options can be added with the `add_option`
    //! function.
    Options() : work_dir_(".") {
        add_option("--work_dir,-w", work_dir_,
                   "Directory containing the MMTF or Hadoop files")
            ->ignore_case()
            ->ignore_underscore()
            ->check(CLI::ExistingDirectory);

        add_option("--ncpu,-n", ncpu_,
                   "Number of CPUs used for run independant jobs")
            ->ignore_case();

        add_option("--entries,-e", entries_,
                   "Preselected index file returned by RCSB")
            ->ignore_case()
            ->check(CLI::ExistingFile);

        add_option("--skip_entries,-s", skip_entries_,
                   "Index file containing PDB ids to skip")
            ->ignore_case()
            ->ignore_underscore()
            ->check(CLI::ExistingFile);
    }

    //! Constructor for an `Options` class which does not use custom options
    //!
    //! This constructor is indended for users who do **not** wish to add custom
    //! options. The *work_dir*, *ncpu*, and *entries* options are added
    //! automatically and the arguments are parsed immediately.
    //! \param argc The number of arguments plus their values plus one.
    //!  Typically obtained from the `main` function.
    //! \param argv The arguments and their values. Typically obtained from the
    //!  `main` function.
    Options(int argc, const char* const argv[]) : Options() { // NOLINT should be C-style
        parse_command_line(argc, argv);
    }

    //! Adds an option to the current workflow
    //!
    //! \param option_name The name of the option on the command line
    //! \param variable The location used to store the variable
    //! \param description The description of the new option
    template<typename T>
    CLI::Option *add_option(std::string option_name, T &variable,
                       std::string description = "") noexcept {
        try {
            return CLI::App::add_option(option_name, variable, description, true);
        } catch (CLI::Error& e) {
            this->exit(e);
            std::exit(1);
        }
    }

    //! Parse the command-line arguments and update containers
    //!
    //! Use this function to read options and update the custom contains if
    //! needed. Once the options are parsed, this function will throw an error.
    //! \param argc The number of arguments plus their values plus one. Typically
    //!  obtained from the `main` function.
    //! \param argv The arguments and their values. Typically obtained from the
    //!  `main` function.
    void parse_command_line(int argc, const char* const argv[]) { // NOLINT should be cy-style
        try {
            parse((argc), (argv));
        } catch (const CLI::Error& e) {
            this->exit(e);
            std::exit(1);
        }
    }

    //! Directory containing the MMTF or Hadoop files
    const std::string& work_dir() const { return work_dir_; }

    //! Number of CPUs used to run independent jobs
    size_t ncpu() const { return ncpu_; }

    //! Index to preselect entries. Eg a search on RCSB
    const std::string& entries() const { return entries_; }

    //! Index to skip entries.
    const std::string& skip_entries() const { return skip_entries_; }

  private:
    std::string work_dir_;
    size_t ncpu_ = 1;
    std::string entries_;
    std::string skip_entries_;
};
} // namespace lemon

#endif
