#ifndef LEMON_HADOOP_HPP
#define LEMON_HADOOP_HPP

#include "lemon/external/gaurd.hpp"

LEMON_EXTERNAL_FILE_PUSH
#ifndef _MSVC_LANG
#include <chemfiles/Frame.hpp>
#include <chemfiles/Trajectory.hpp>

#include <arpa/inet.h>
#include <dirent.h>
#else
#include <WinSock2.h>
#include <chemfiles/Frame.hpp>
#include <chemfiles/Trajectory.hpp>
#include <lemon/external/dirent.hpp>
#endif
LEMON_EXTERNAL_FILE_POP

#include <cassert>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include <array>

namespace lemon {

//! The `Hadoop` class is used to read input sequence files.
//!
//! This class reads an Apache Hadoop Sequence file and iterates through the
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
    //! Be-warned that this function does minimal error checking and should only
    //! be used if has_next() has returned `true`.
    //! \return A pair of `std::vector<char>`s. The first member contains the
    //! PDB ID and the second contains the GZ compressed MMTF file.
    std::pair<std::string, std::vector<char>> next() { return read(); }

    //! The size of the starting header
    static auto constexpr HADOOP_HEADER_SIZE = 90;
  private:
    std::istream& stream_;
    std::string marker_ = "";
    std::vector<char> key_;

    // Initialize the sequence file.
    void initialize_() {
        // Completly skip the header as it is the same in all RCSB Hadoop files
        stream_.exceptions(std::ifstream::badbit | std::ifstream::failbit);
        std::array<char, HADOOP_HEADER_SIZE> buffer;
        stream_.read(buffer.data(), HADOOP_HEADER_SIZE - 3);
    }

    // Read four bytes and return as an 4 byte integer
    int read_int() {
        int ret;
        stream_.read(reinterpret_cast<char*>(&ret), 4);
        return static_cast<int>(ntohl(static_cast<uint32_t>(ret)));
    }

    std::pair<std::string, std::vector<char>> read() {
        auto sync_check = read_int();
        auto constexpr MARKER_SIZE = 16;

        if (sync_check == -1) {
            std::vector<char> marker(MARKER_SIZE);
            stream_.read(marker.data(), MARKER_SIZE);
            // Only valid if using the full version
            // assert(std::string(marker.data(), 16) == marker_);
            return this->read();
        }

        auto key_length = read_int();

        // Do not check this during runtime as it should all be the same
        assert(key_length >= 4);

        std::vector<char> key(static_cast<size_t>(key_length));
        stream_.read(key.data(), key_length);
        const auto entry = std::string(key.data() + 1, 4);

        assert(sync_check >= 8);

        // Remove junk characters added by Java Serialization
        std::array<char, 4> junk;
        stream_.read(junk.data(), 4);

        int value_length = sync_check - key_length;
        std::vector<char> value(static_cast<size_t>(value_length - 4));
        stream_.read(value.data(), value_length - 4);

        return {entry, value};
    }
};

//! \brief Read a directory containing hadoop sequence files
inline std::vector<std::string> read_hadoop_dir(const std::string& p) {

    DIR* dp;
    std::vector<std::string> pathvec;
    pathvec.reserve(700); // NOLINT typlically 700. Keeping the magic number

    dp = opendir(p.c_str());
    if (dp == nullptr) {
        throw std::runtime_error("Path does not exist or could not be read.");
    }

    for (auto entry = readdir(dp); entry != nullptr; entry = readdir(dp)) {
        if (entry->d_name[0] == '_' || entry->d_name[0] == '.' ||
            entry->d_type != DT_REG) {
            continue;
        }

        std::string s = entry->d_name;
        if (s.find('.') != std::string::npos) {
            closedir(dp);
            throw std::runtime_error(
                "Directory provided has file with extensions.\nPlease "
                "remove files with extensions if you are sure the seqeunce "
                "files are valid.");
        }
        s.insert(0, "/");
        s.insert(0, p);
        pathvec.emplace_back(s);
    }

    closedir(dp);
    return pathvec;
}
} // namespace lemon

#endif
