#ifndef LEMON_OPTIONS_HPP
#define LEMON_OPTIONS_HPP

#include <boost/program_options.hpp>

namespace lemon {

namespace po = boost::program_options;

//! The `Options` class is used to read command line arguments.
//!
//! A majority of **Lemon** workflows begin with the parsing of various command-
//! line arguments. This class helps the user parse the most common arguments
//! such as the directory containing the Hadoop sequence files, the number of
//! cores to use for a given search or a preselected set of entries.
//! Users can add their own options with the `add_options` member function.
class Options {
   public:
    //! Default constructor to initialize a custom `Options` class
    //!
    //! This constructor is intended for users who wish to add their own custom
    //! options. The *work_dir*, *ncpu*, and *entries* options are added
    //! automatically and additional options can be added with the `add_option`
    //! function.
    Options()
        : desc_("Options"), ncpu_(1), work_dir_("."), options_compiled_(false) {
        add_option("work_dir,w", work_dir_,
                   "Directory containing the MMTF or Hadoop files");

        add_option("ncpu,n", ncpu_,
                   "Number of CPUs used for run independant jobs");

        add_option("entries,e", entries_,
                   "Preselected index file returned by RCSB");
    }

    //! Constructor for an `Options` class which does not use custom options
    //!
    //! This constructor is indended for users who do **not** wish to add cusom
    //! options. The *work_dir*, *ncpu*, and *entries* options are added
    //! automatically and the arguments are parsed immediately.
    //! \param argc The number of argments plus their values plus one. Typically
    //!  obtained from the `main` function.
    //! \param argv The arguments and their values. Typically obtained from the
    //!  `main` function.
    Options(int argc, const char* const argv[]) : Options() {
        parse_command_line(argc, argv);
    }

    //! Add a custom option
    //!
    //! Use this member function to add additional options to the current
    //! `Options` class. The `container` argument is updated when
    //! `parse_command_line` is called. Once the options are parsed, this
    //! function will throw an error as it is not possible to add arguments
    //! after parsing. \param name The name of argument. Use a comma to add a
    //! shortcut. For
    //!   example, `ncpu,n` will recognize `--ncpu` and `-n`.
    //! \param container The variable where the command's argument will be
    //! stored. \param description Optional additional help for the user.
    template <typename T>
    void add_option(const std::string& name, T& container,
                    const std::string& description = "") {
        if (options_compiled_) {
            throw std::runtime_error("You cannot add options after parsing.");
        }
        desc_.add_options()(name.c_str(),
                            po::value<T>(&container)->default_value(container),
                            description.c_str());
    }

    //! Parse the command-line arguments and update containers
    //!
    //! Use this function to read options and update the custom contains if
    //! needed. Once the options are parsed, this function will throw an error.
    //! \param argc The number of argments plus their values plus one. Typically
    //!  obtained from the `main` function.
    //! \param argv The arguments and their values. Typically obtained from the
    //!  `main` function.
    void parse_command_line(int argc, const char* const argv[]) {
        if (options_compiled_) {
            throw std::runtime_error("You cannot parse options twice.");
        }
        po::store(po::parse_command_line(argc, argv, desc_), vm_);
        po::notify(vm_);
        options_compiled_ = true;
    }

    //! Directory containing the MMTF or Hadoop files
    const std::string& work_dir() { return work_dir_; }

    //! Number of CPUs used for run independant jobs
    size_t ncpu() const { return ncpu_; }

    //! Index to preselect entries. Eg a search on RCSB
    const std::string& entries() { return entries_; }

   private:
    po::options_description desc_;
    po::variables_map vm_;

    bool options_compiled_;
    std::string work_dir_;
    size_t ncpu_;
    std::string entries_;
};
}  // namespace lemon

#endif
