#ifndef HADOOP_HPP
#define HADOOP_HPP

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
#include <future>

namespace lemon {

namespace fs = boost::filesystem;

class Hadoop {
   public:
    Hadoop(std::istream& stream) : stream_(stream) { initialize_(); }

    bool has_next() { return stream_.peek() != std::char_traits<char>::eof(); }

    std::pair<std::vector<char>, std::vector<char>> next() { return read(); }

   private:
    std::istream& stream_;
    std::string marker_ = "";
    std::vector<char> key_;

    void initialize_() {
        // Completly skip the header as it is the same in all RCSB Hadoop files
        stream_.exceptions(std::ifstream::badbit | std::ifstream::failbit);
        char buffer[86];
        stream_.read(buffer, 87);
    }

    int read_int() {
        int ret;
        stream_.read(reinterpret_cast<char*>(&ret), 4);
        return ntohl(ret);
    }

    std::pair<std::vector<char>, std::vector<char>> read() {
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

        assert(sync_check >= 8);

        // Remove junk characters added by Java Serialization
        char junk[4];
        stream_.read(junk, 4);

        int value_length = sync_check - key_length;
        std::vector<char> value(value_length - 4);
        stream_.read(value.data(), value_length - 4);

        return {key, value};
    }
};

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
