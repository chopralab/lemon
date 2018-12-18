#ifndef LEMON_HADOOP_HPP
#define LEMON_HADOOP_HPP

#include <chemfiles/Frame.hpp>
#include <chemfiles/Trajectory.hpp>

#ifndef _MSVC_LANG
#include <arpa/inet.h>
#else
#include <WinSock2.h>
#endif

#include <cassert>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

namespace lemon {

namespace fs = boost::filesystem;

//! The `Hadoop` class is used to read input sequence files.
//!
//! This class reads an Apache Hadoop Sequence file and interates through the
//! key/value pairs.  It has been modified so that it can only read files
//! supplied by RCSB at this location:
//! [Full](https://mmtf.rcsb.org/v1.0/hadoopfiles/full.tar)
class Hadoop {
   public:

    //! Create a `Hadoop` class using a `std::istream`.
    //!
    //! This Hadoop constructor takes a binary stream as input.  This stream
    //! must be open and contain data from a sequence file obtained from RCSB.
    Hadoop(std::istream& stream) : stream_(stream) { initialize_(); }

    //! Returns if a sequence file has remaining MMTF records in it.
    //!
    //! Use this function to check if the sequence file has any remaining MMTF
    //! records stored in it.
    //! \return True if another MMTF record is present. False otherwise.
    bool has_next() { return stream_.peek() != std::char_traits<char>::eof(); }

    //! Returns the next MMTF file.
    //!
    //! This function reads the next MMTF record from the underlying stream.
    //! Bewarned that this function does minimal error checking and should only
    //! be used if has_next() has returned `true`.
    //! \return A pair of `std::vector<char>`s. The first member contains the
    //! PDB ID and the second contains the GZ compressed MMTF file.
    std::pair<std::string, std::vector<char>> next() { return read(); }

   private:
    std::istream& stream_;
    std::string marker_ = "";
    std::vector<char> key_;

    // Initialize the sequence file.
    void initialize_() {
        // Completly skip the header as it is the same in all RCSB Hadoop files
        stream_.exceptions(std::ifstream::badbit | std::ifstream::failbit);
        char buffer[90];
        stream_.read(buffer, 87);
    }

    // Read four bytes and return as an 4 byte integer
    int read_int() {
        int ret;
        stream_.read(reinterpret_cast<char*>(&ret), 4);
        return ntohl(ret);
    }

    std::pair<std::string, std::vector<char>> read() {
        auto sync_check = read_int();

        if (sync_check == -1) {
            std::vector<char> marker(16);
            stream_.read(marker.data(), 16);
            // Only valid if using the full version
            // assert(std::string(marker.data(), 16) == marker_);
            return this->read();
        }

        auto key_length = read_int();

        // Do not check this during runtime as it should all be the same
        assert(key_length >= 4);

        std::vector<char> key(key_length);
        stream_.read(key.data(), key_length);
        const auto entry = std::string(key.data() + 1, 4);

        assert(sync_check >= 8);

        // Remove junk characters added by Java Serialization
        char junk[4];
        stream_.read(junk, 4);

        int value_length = sync_check - key_length;
        std::vector<char> value(value_length - 4);
        stream_.read(value.data(), value_length - 4);

        return {entry, value};
    }
};

//! \brief Read a directory containing hadoop sequence files
inline std::vector<fs::path> read_hadoop_dir(const fs::path& p) {
    if (!fs::is_directory(p)) {
        throw std::runtime_error("Provided directory not valid.");
    }

    std::vector<fs::path> pathvec;
    pathvec.reserve(600);

    // There's only ~550 files to read here!
    auto begin = fs::directory_iterator(p);
    fs::directory_iterator end;
    std::transform(
        begin, end, std::back_inserter(pathvec),
        [](fs::directory_entry& entry) {
            if (fs::is_directory(entry.path())) {
                throw std::runtime_error(
                    "Directory provided has subdirectories.\nPlease make sure "
                    "you are using the tar ball provided by RCSB.");
            }
            if (entry.path().has_extension()) {
                throw std::runtime_error(
                    "Directory provided has file with extensions.\nPlease "
                    "remove files with extensions if you are sure the seqeunce "
                    "files are valid.");
            }
            return entry.path();
        });

    return pathvec;
}
}  // namespace lemon

#endif
