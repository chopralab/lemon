#ifndef LEMON_OPTIONS_HPP
#define LEMON_OPTIONS_HPP

#include <boost/program_options.hpp>

namespace lemon {

namespace po = boost::program_options;

/*!
 * The `Options` class is used to read command line arguments.
*/
class Options {
   public:
    Options(int argc, const char* const argv[]) {
        po::options_description desc("Options");
        desc.add_options()(
            "work_dir,w", po::value<std::string>()->default_value("."),
            "Directory containing the MMTF or Hadoop files")(
            "distance,d", po::value<double>()->default_value(6.0),
            "Distance cutoff used for distance based searches")(
            "reference,r", po::value<std::string>()->default_value(""),
            "Reference file for structure based searches")(
            "ncpu,n", po::value<int>()->default_value(1),
            "Number of CPUs used for run independant jobs")(
            "entries,e", po::value<std::string>()->default_value(""),
            "Index file returned by RCSB used");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        work_dir_ = vm["work_dir"].as<std::string>();
        distance_ = vm["distance"].as<double>();
        reference_ = vm["reference"].as<std::string>();
        ncpu_ = vm["ncpu"].as<int>();
        entries_ = vm["entries"].as<std::string>();
    }

    //! Directory containing the MMTF or Hadoop files
    const std::string& work_dir() {
        return work_dir_;
    }

    //! Distance cutoff used for distance based searches
    double distance() const {
        return distance_;
    }

    //! Reference file for structure based searches
    const std::string& reference() const {
        return reference_;
    }

    //! Number of CPUs used for run independant jobs
    size_t ncpu() const {
        return ncpu_;
    }

    //! Index file returned by RCSB
    const std::string& entries() {
        return entries_;
    }

   private:

    std::string work_dir_;

    double distance_;

    std::string reference_;

    size_t ncpu_;

    std::string entries_;
};
}

#endif
